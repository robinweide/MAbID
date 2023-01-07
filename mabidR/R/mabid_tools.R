#' mabid_bind
#'
#' `mabid_bind()` binds TTAA and GATC SCE's
#'
#' This is the single-cell loader of AbID-data. It will do normalization on control and
#' mappability and generate multiple handy tracks and metadata.
#'
#' @param sce_1 A SCE object
#' @param sce_2 A SCE object
#' @return An SCE-object of the RangedSummarizedExperiment-subclass (logcounts = TMM-normalised lfc)
#' @export
mabid_bind  <- function(sce_1, sce_2){
  # check colnames
  cn <- c(colnames(sce_1), colnames(sce_2))
  if(length(cn) != length(unique(cn))){
    colnames(sce_1) <- paste0(colnames(sce_1),".x")
    colnames(sce_2) <- paste0(colnames(sce_2),".y")
  }

  # get names
  RN = table(c(rownames(sce_1),rownames(sce_2)))
  RN <- names(RN[RN == 2])

  CN = table(c(colnames(colData(sce_1)),
               colnames(colData(sce_2))))
  CN <- names(CN[CN == 2])

  # massage rowdata
  sce_1_rd <- rowData(sce_1[RN,])
  sce_2_rd <- rowData(sce_2[RN,])

  RD <- DataFrame(control = rowSums(cbind(sce_1_rd$control/sum(sce_1_rd$control, na.rm = T),
                                          sce_2_rd$control/sum(sce_2_rd$control, na.rm = T)))*1e6,
                  mappability = rowSums(cbind(sce_1_rd$mappability/sum(sce_1_rd$mappability, na.rm = T),
                                              sce_2_rd$mappability/sum(sce_2_rd$mappability, na.rm = T)))*1e6,
                  row.names = rownames(sce_1_rd))

  sce_1 <- sce_1[RN, ]
  sce_2 <- sce_2[RN, ]
  rowData(sce_1) <- RD
  rowData(sce_2) <- RD

  # massage coldata
  colData(sce_1) <- colData(sce_1)[,CN]
  colData(sce_2) <- colData(sce_2)[,CN]

  out = cbind(sce_1,
              sce_2)
  out$ID <- colnames(out)
  out
}

#' mabid_rebin
#'
#' `mabid_rebin()` re-bins a SCE
#'
#' Change the bin-size a sce-object containing MAbID-datasets.
#'
#' @param sce Object: SCE
#' @param bin_size Integer: size of a bin in bp
#' @param threads Integer: number of CPU's
#' @return An SCE-object of the RangedSummarizedExperiment-subclass (logcounts = TMM-normalised lfc)
#' @export
mabid_rebin <- function(sce, bin_size = 20e3, threads = 7){
  pseudocount <- sce$pseudocount[1]
  splitted <- split(sce ,seqnames(sce))
  per_chrom <- mclapply(mc.cores = threads, splitted, function(scec){
    idx <-   floor(BiocGenerics::start(scec) / bin_size) * bin_size

    new_counts   <- aggregate(as.matrix(BiocGenerics::counts(scec)),by=list(idx), sum)
    new_controls <- aggregate(SummarizedExperiment::rowData(scec)$control,by=list(idx), sum)
    new_mapp     <- aggregate(SummarizedExperiment::rowData(scec)$mappability,by=list(idx), sum)

    mat_control <- new_controls[,-1]
    mat_counts <- as.matrix(new_counts[,-1])
    mat_mappability  <- new_mapp[,-1]

    if(!'matrix' %in% class(mat_counts)){
      mat_counts <- matrix(mat_counts,ncol = 1)
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

    y <- SingleCellExperiment::SingleCellExperiment(list(counts = Matrix::Matrix(mat_counts,sparse = T),
                                                         logcounts = Matrix::Matrix(mat_oec,sparse = T),
                                                         tpm = Matrix::Matrix(RPM_signal,sparse = T)
    ))
    SummarizedExperiment::rowData(y)$control     <- mat_control
    SummarizedExperiment::rowData(y)$mappability <- mat_mappability
    SummarizedExperiment::rowRanges(y) <- GenomicRanges::GRanges(as.character(seqnames(scec)[1]),
                                                                 IRanges::IRanges(start = new_counts[,1],
                                                                                  width = bin_size))
    rownames(y) <- as.character(SummarizedExperiment::rowRanges(y))
    y
  })

  out <- suppressMessages(do.call('rbind',per_chrom))

  SummarizedExperiment::colData(out)  <- SummarizedExperiment::colData(sce)[colnames(out),]
  out$bin_size <- bin_size

  GenomeInfoDb::seqinfo(out) <- GenomeInfoDb::seqinfo(sce)
  # done!
  out

}


#' mabid_IC
#'
#' `mabid_IC()` calculates and adds information content
#'
#'
#' @param sce A SCE object
#' @param counts_mappability A named vector of mappability-counts. Names should be in 'chrX:0-1' format
#' @param p Number of CPU's to use
#' @return An SCE-object of the RangedSummarizedExperiment-subclass, with an added "IC" colData-column
#' @export
mabid_IC <- function(sce, counts_mappability, p = 7){

  counts_mappability <- counts_mappability[match(rownames(sce), names(counts_mappability))]
  cells <- 1:ncol(sce)
  cells <- cells[sce$depth > 0]
  out <- rbindlist(mclapply(mc.cores = p, X = cells, function(i){  # loop over samples
    Xi1 <- as.vector(counts(sce[,i]))
    Xij = cbind(Xi1,counts_mappability)

    # smooth
    Xij <- apply(Xij, 2, function(e){ksmooth(1:length(e),e, kernel = 'normal',n.points = length(e),bandwidth = 6.334)$y})

    # rm 0 bins
    Xij <- Xij[Xij[,1] > 0,]
    Xij_no0 <- Xij
    # d is the total sum of Xij
    d <- sum(Xij)

    # Pi is the marginal probability of observing counts in bin i
    Pi <- apply(Xij,1, function(e){ sum(e)/d })

    # Qj is the marginal probability of observing counts in sample j
    Qj <- apply(Xij,2, function(e){ sum(e)/d })

    # E is the expected number of counts in bin i in sample j
    Eij <- (d*cbind(Pi * Qj[1],Pi * Qj[2]))

    # GoF of our predictions compared to the actual counts per bin:
    g <- rowSums(cbind((Xij[,1] - Eij[,1])/Eij[,1],
                       (Xij[,2] - Eij[,2])/Eij[,2]))

    # iterate to convergence
    GOF_range <- quantile(g,c(0.05,0.5))
    S <- rep(T,length(g))
    Xij_sub <- NA
    for(e in 1:20){
      S_update = (g >= GOF_range[1]) & (g < GOF_range[2])

      if(sum(S_update != S) < 5){ break }

      S <- S_update

      Xij_sub <- Xij[S_update,]
      # d is the total sum of Xij
      d <- sum(Xij_sub)

      # Qj is the marginal probability of observing counts in sample j
      Qj <- apply(Xij_sub, 2, function(e){ sum(e)/d })

      # E is the expected number of counts in bin i in sample j
      E <- (d*cbind(Pi * Qj[1],Pi * Qj[2]))

      # GoF of our predictions compared to the actual counts per bin:
      g <- rowSums(cbind((Xij[,1] - E[,1])/E[,1],
                         (Xij[,2] - E[,2])/E[,2]))

      GOF_range <- quantile(g,c(0.05,0.5))
    }

    a = colSums(Xij_sub)/sum(Xij_sub) # final
    b = colSums(Xij_no0)/sum(Xij_no0) # initial

    data.table('ID' = sce$ID[i], 'IC' = (a[2]/b[2])/(a[1]/b[1]))
  }))

  # add entries wher depth was zero to IC0
  out <- rbind(out,data.table('ID' = sce$ID[!sce$ID %in% out$ID], 'IC' = 0))

  # add IC-column to sce
  sce <- mabid_annotate(sce, anno_vec = setNames(out$IC,out$ID), where = 'ID', name = 'IC')

  # return sce
  sce
}
