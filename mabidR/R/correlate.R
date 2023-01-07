#
#' General correlation function
#'
#' `mabid_correlate()`
#'
#' Get different correlation-metrics (region and chromosome-wide)
#'
#' @param sce A SCE object
#' @param blacklist Ignore these regions
#' @param ROI Region of interest to correlate to
#' @param GRlist_other add a list of bigwig-GRs to correlate to
#' @param assay_to_use assay in sce to use
#' @param bin_size bin-size to use
#' @param kernel_size Smoothing window
#' @param winsorize Winsorize (Replace Extreme Values by Less Extreme Ones)
#' @param npc threads to use
#' @param ignore_negatives set negative values to zero
#' @export
mabid_correlate <- function(sce, blacklist, ROI = NULL, GRlist_other = NULL,
                            assay_to_use = 'logcounts', bin_size = NULL,
                            kernel_size = 4, winsorize = 0.9995, npc = 7, ignore_negatives = F){
  if(is.null(bin_size)){
    bin_size <- median(sce$bin_size)
  }

  #### REGION
  if(is.null(ROI)){
    ROI <- seqlengths(sce)
    ROI <- GRanges(names(ROI), IRanges(1,ROI))
    seqinfo(ROI) <- seqinfo(sce)
  }


  #### MAbID
  dat <- assay(subsetByOverlaps(subsetByOverlaps(sce,ROI), blacklist, invert = T),assay_to_use)
  dat <- data.table(reshape2::melt(as.matrix(dat)))
  dat <- dat[,  c('seqnames', 'start', 'end') := tstrsplit(Var1, "[-:]")][]
  dat$start <- as.numeric(dat$start)
  dat$end <- as.numeric(dat$end)

  dat <- rbindlist(lapply(split(dat, list(dat$Var2, dat$seqnames)), function(dt){
    dt$value <- DescTools::Winsorize(dt$value, probs = c(1-winsorize,winsorize))
    dt
  }))

  colnames(dat)[1:2] <- c('binID','sample')
  dat$binID <- NULL

  #### other
  if(!is.null(GRlist_other)){
    dat$sample <- paste0(dat$sample,'_MAbID')

    datO <- data.table(as.data.frame(rbindlist(lapply(GRlist_other, function(e){
      e <- as.data.frame(subsetByOverlaps(subsetByOverlaps(e,ROI), blacklist, invert = T))
      e$value <- DescTools::Winsorize(e$score, probs = c(1-winsorize,winsorize))
      e
    }),idcol = 'sample')))
    datO$sample <- paste0(datO$sample,'_other')

    dat <- rbind(dat, datO[,colnames(dat), with = F])
  }

  #### cor-tests
  dat$pos <- apply(cbind(dat$seqnames,floor(dat$start/bin_size)*bin_size),
                   1,
                   paste0, collapse = ":")

  kill_inf <- function(x){
    x[x ==  Inf] <- max(x[x !=  Inf])
    x[x ==  -Inf] <- min(x[x != -Inf])
    x
  }
  dat <- dat[,list(value = kill_inf(value)), by = 'sample,pos,seqnames,start']

  if(ignore_negatives){
    dat$value[dat$value < 0] <- 0
  }
  # fig1c-like
  mat <- dat[,mean(value), by = 'sample,pos']
  mat <- as.data.frame(dcast(mat, pos ~ sample, value.var = 'V1'))
  rownames(mat) <- mat$pos; mat$pos <- NULL

  f1 <- apply(mat, 2, scale)
  f1_p <- reshape2::melt(cor(f1,method = 'p',use = "pairwise.complete.obs"))
  f1_s <- reshape2::melt(cor(f1,method = 's',use = "pairwise.complete.obs"))
  out_f1 <- data.table(merge(f1_p,f1_s, by = 1:2))
  colnames(out_f1) <- c('A','B','R_fig1', 'rho_fig1')

  # fig4-like
  mat <- rbindlist(lapply(split(dat, list(dat$sample, dat$seqnames)), function(dt){
    dt <- dt[order(dt$start)]
    dt$value <- ksmooth(dt$start, dt$value, kernel = "normal",
                        bandwidth = kernel_size*bin_size, n.points = nrow(dt))$y
    dt
  }))
  mat <- mat[,mean(value), by = 'sample,pos']
  mat <- as.data.frame(dcast(mat, pos ~ sample, value.var = 'V1'))
  rownames(mat) <- mat$pos; mat$pos <- NULL

  f4_p <- reshape2::melt(cor(mat,method = 'p',use = "pairwise.complete.obs"))
  f4_s <- reshape2::melt(cor(mat,method = 's',use = "pairwise.complete.obs"))
  out_f4 <- data.table(merge(f4_p,f4_s, by = 1:2))
  colnames(out_f4) <- c('A','B','R_fig4', 'rho_fig4')

  out <- merge(out_f1,out_f4, by = c('A','B'))

  if(ignore_negatives){
    colnames(out)[3:6] <- paste0(colnames(out)[3:6],'_pos')
  }

  return(out)
}
