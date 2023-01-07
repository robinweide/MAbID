get_OE <- function(signal, control, mapp0_2_oe0 = T){
  OE_links <- (signal + 1) / (control + 1)

  ToBw <- sum(signal, na.rm = T) + (length(signal) + 1)
  TeBw <- sum(control, na.rm = T) + (length(control) + 1)

  OE_rechts <- TeBw/ToBw

  OE <- OE_links * OE_rechts

  # set OE to 1 if signal == 0 AND control == 0
  OE[signal == 0 & control == 0] <- 1

  if(mapp0_2_oe0){
    OE[control == 0] <- 0
  }

  OE
}
#' Load an single-cell MAb-ID-sample
#'
#' `load_mabid()` returns a SingleCellExperiment-object.
#'
#' This is the single-cell loader of AbID-data. It will do normalization on control and
#' mappability and generate multiple handy tracks and metadata.
#'
#' @param signal_hdf5 Path(s) to a hdf5-file (binned or per TTAA).
#' @param map_hdf5 Path to a matching mappability.hdf5
#' @param control_hdf5 Path(s) to matching control.hdf5 (typically "nones")
#' @param reference A BSgenome-object
#' @param bin_size The used bin size of the signal-hdf5
#' @param sample_name String to for colnaames(sce). Go for a discriptive one!
#' @param threads number of threads to use
#' @param blacklist A GRanges-object with regions that should be set to NA due to them being not representative. Also include (peri)centromeric regions!
#' @param population_average make a population_average?
#' @param pseudocount Default = 1
#' @param min_bin skip small chromosomes (usually patches)
#' @param downsample Downsample counts of all samples to lowest (can also be a number)
#' @param downsample_control Downsample control to a number
#' @return An SCE-object of the RangedSummarizedExperiment-subclass (logcounts = TMM-normalised lfc, normcounts = GLM-residuals)
#' @examples
#' \dontrun{
#' output_sce <- load_scMAbid(signal_hdf5 = "index_30.ABBC_001.SBC_001.100000.hdf5",
#'                        control_hdf5 = "index_30.ABBC_001.SBC_007.100000.hdf5",
#'                        map_hdf5 = 'Homo_sapiens.GRCh37.dna.primary_assembly.TTAA.100000.hdf5',
#'                        bin_size = 100000, reference = BSgenome.Hsapiens.UCSC.hg19)
#' }
#' @export
load_mabid <- function(signal_hdf5, map_hdf5, control_hdf5, reference, bin_size = NULL,
                       sample_name = NULL, threads = 5, blacklist = NULL,
                       population_average = F, pseudocount = 1, min_bin = 100,
                       downsample = F,downsample_control = NULL){

  # changelog:
  ## 28/07/22: EpiDamID RPKM-normalisation and log2OE method
  ## 16/08/22: added downsampling

  require(SingleCellExperiment)

  if(is.null(sample_name)){
    sample_name <- gsub('.hdf5$','',BiocGenerics::basename(signal_hdf5))
  }

  signal_hdf5 <- c(signal_hdf5, map_hdf5, control_hdf5)
  H5LS <- rhdf5::h5ls(signal_hdf5[1])
  H5LS <- H5LS[as.numeric(H5LS$dim) >= min_bin,]

  GR <- bin_size*as.numeric(H5LS$dim)
  names(GR) <- H5LS$name

  GR = GenomicRanges::tileGenome(GR,
                                 tilewidth = bin_size, cut.last.tile.in.chrom = T)

  # convert to the same scheme as the reference uses
  newStyle <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(GR), GenomeInfoDb::provider(reference))
  GR  <- GenomeInfoDb::renameSeqlevels(GR, na.omit(newStyle))
  chroms_to_use <- GenomeInfoDb::seqlevels(GR)[GenomeInfoDb::seqlevels(GR) %in% GenomeInfoDb::seqlevels(reference)]
  GR <- GenomeInfoDb::keepSeqlevels(GR, chroms_to_use, pruning.mode = 'coarse')

  rr = GenomicRanges::GRanges(GenomeInfoDb::seqinfo(reference))
  rr <- GenomeInfoDb::keepSeqlevels(rr, chroms_to_use, pruning.mode = 'coarse')

  suppressWarnings(GenomeInfoDb::seqinfo(GR) <- GenomeInfoDb::seqinfo(rr))
  GR <- GenomicRanges::trim(GR)


  chroms_to_use <- H5LS$name[H5LS$name %in% names(newStyle[!is.na(newStyle)])]

  # get counts
  cnts <- parallel::mclapply(mc.cores = threads, chroms_to_use, function(chrom){
    mat <- lapply(signal_hdf5, rhdf5::h5read, name = chrom)
    mat <- sapply(mat, function(x){ x[1:length(mat[[1]])]})
    mat <- as.matrix(mat)
    colnames(mat) <- c(sample_name, 'mappability','control')
    mat
  })
  cnts <- do.call('rbind',cnts)

  mat_counts      <- cnts[,-c(ncol(cnts)-1,ncol(cnts))]
  mat_control     <- cnts[,c(ncol(cnts))]
  mat_mappability <- cnts[,c(ncol(cnts)-1)]

  if(!is.null(blacklist)){
    if(class(blacklist) == 'GRanges'){
      idx <- unique(GenomicRanges::findOverlaps(GR,blacklist)@from)
      mat_mappability[idx] <- 0
    } else {
      stop('blacklist should be a GRanges-object!')
    }
  }

  mat_counts[mat_mappability == 0] <- 0
  mat_control[mat_mappability == 0] <- 0

  if(!'matrix' %in% class(mat_counts)){
    mat_counts <- matrix(mat_counts,ncol = 1)
  }

  if(population_average){
    mat_counts = matrix(rowSums(mat_counts,na.rm = T), ncol = 1)
    sample_name <- 'population_average'
  }

  if(!downsample == F){
    target <- min(colSums(mat_counts))
    if(is.numeric(downsample)){
      target <- downsample
    }
    mat_counts <- as.matrix(Seurat::SampleUMI(mat_counts, max.umi = target))
  }

  if(!is.null(downsample_control)){
    mat_control <- as.vector(Seurat::SampleUMI(Matrix::Matrix(mat_control,ncol = 1,sparse = T), max.umi = downsample_control)[,1])
  }

  # RPKM Log2(O/E) normalise
  SF_control <- sum(mat_control)/1e6
  SF_signal  <- colSums(mat_counts)/1e6

  RPM_control <- mat_control/SF_control
  RPM_signal  <- t(t(mat_counts)/SF_signal)

  SFb <- (bin_size/1e3)

  RPKM_control <- RPM_control / SFb
  RPKM_signal  <- RPM_signal / SFb

  RPKM_control <- RPKM_control + pseudocount
  RPKM_signal   <- RPKM_signal + pseudocount

  mat_oec = log2(RPKM_signal/RPKM_control)
  mat_oec[mat_mappability == 0] <- 0

  # get basic stats for coldata
  sp = signal_hdf5[-c(ncol(cnts)-1,ncol(cnts))]
  if(population_average){
    sp = paste0(sp, collapse = ",")
  }
  cdata <- S4Vectors::DataFrame(signal_path      = sp,
                                control_path     = control_hdf5,
                                mappability_path = map_hdf5,
                                bin_size         = bin_size,
                                depth            = colSums(mat_counts,na.rm = T),
                                pseudocount      = pseudocount,
                                reference        = reference@pkgname)

  # make sce
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = Matrix::Matrix(mat_counts,sparse = T),
                                                         logcounts = Matrix::Matrix(mat_oec,sparse = T),
                                                         tpm = Matrix::Matrix(RPM_signal,sparse = T)
  ))
  SummarizedExperiment::colData(sce) <- cdata
  colnames(sce) <- sample_name

  SummarizedExperiment::rowRanges(sce) <- GR
  rownames(sce) <- as.character(GR)

  rowData(sce)$control     <- mat_control
  rowData(sce)$mappability <- mat_mappability

  # done!
  sce
}
