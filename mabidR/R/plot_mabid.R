#' Pull values of a ROI from the SCE.
#'
#' `mabid_get_track()` returns a datatable of the values of a regions of interest
#'
#' This enables quick querying and subsetting of an SCE
#'
#' @param sce An MAb-ID SCE object
#' @param chrom Chromosome of interest
#' @param start Start of region
#' @param end End of region
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param metadata Will add all colData. Can slow things down!!!
#' @return A data.table-object
#' @examples
#' \dontrun{
#' data_roi <- m_get_track(sce_SBC001,
#'                            chrom = 'chr1', start = 0, end = 100e6,
#'                            assay_to_use = 'OE_control')
#' }
#' @export
mabid_get_track <- function(sce, chrom = 'chr1', start = 0, end = 5000000, assay_to_use = 'logcounts', metadata = T){
  ROI <- GenomicRanges::GRanges(paste0(chrom,":",start,'-',end,":*"))
  mini_sce <- IRanges::subsetByOverlaps(sce, ROI)
  dt <- SummarizedExperiment::assay(mini_sce,assay_to_use)
  vals <- lapply(1:ncol(dt), function(i){dt[,i]})
  names(vals) <- colnames(dt)
  vals <- data.table::rbindlist(lapply(1:length(vals), function(i){
    data.table::data.table( x = names(vals[[i]]), sample = names(vals[i]), value = vals[[i]])
  }))
  dt <- cbind(reshape2::colsplit(vals$x,pattern = '[:-]',c('seqnames','start','end')), data.table::data.table(vals[,-1]))
  dt$variable <- assay_to_use

  if(metadata){
    mt = colData(sce)
    mt$sample = rownames(mt)
    dt = merge(dt, mt, by = 'sample')
  }
  dt
}

#' Plot values of a ROI from a SCE as barplots.
#'
#' `mabid_trackplot()` returns generates barplots, coloured on sample.
#'
#' Pretty pretty!
#'
#' @param sce An MAb-ID SCE object
#' @param chrom Chromosome of interest
#' @param start Start of region
#' @param end End of region
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param GR add bigwig-like tracks.
#' @param colour_on Use a column of the output of `mabid_get_track(sce)` (default: 'sample')
#' @param kernelLen (integer) Number of data-points to consider in the kernel.
#' @param ylim Set range of y-axis
#' @param ylim_per_sample Automagically get nice ylims
#' @param annot_GRlist a named list of GR's
#' @param annot_col Set a standard colour for de annot_GRlist
#' @param hline Add a horizontal line at a y-axis value
#' @param plot_ratio Set aspect ratio (default: 1/10)
#' @return A Ggplot2-object
#' @examples
#' \dontrun{
#' data_roi <- mabid_get_track(sce_SBC001,
#'                            chrom = 'chr1', start = 0, end = 100e6,
#'                            assay_to_use = 'OE_control')
#' }
#' @export
mabid_trackplot <- function(sce, chrom = 'chr1', start = 0, end = 5000000, colour_on = 'sample',kernelLen = 3,
                            annot_GRlist = NULL,annot_col = NULL,
                            assay_to_use = 'logcounts', GR = NULL,
                            ylim = NULL, ylim_per_sample = T, hline = NULL,plot_ratio = 1/10){
  floor_dec <- function(x, level=2) round(x - 5*10^(-level-1), level)

  legend_label <- switch(assay_to_use,
                         'OE_control_log2' = expression(Log[2]~(frac(counts, control))),
                         'OE_mappability_log2' = expression(Log[2]~(frac(counts, mappability))),
                         'counts' = 'counts','normcounts' = 'normalised signal',
                         'logcounts' = expression(Log[2]~(frac(counts, control))),
                         assay_to_use)

  if(!is.null(annot_GRlist)){
    annot <- data.table::rbindlist(lapply(annot_GRlist,function(x){
      as.data.frame(x[seqnames(x) == chrom])
    }),idcol = 'sample',)
    annot <- as.data.frame(annot)
    annot <- annot[annot$end >= start & annot$start <= end,]

    annot$start[annot$start < start] <- start
    annot$end[annot$end > end] <- end
  }

  dat <- mabid_get_track(sce,
                         chrom = chrom,start = start,end = end,assay_to_use = assay_to_use)
  dat_tmp <- dat
  dat <- tidyr::spread(dat[,c('sample','start','value')],key = 'sample',value = 'value')
  dat <- dat[order(dat$start),]
  samplingTime = median(sce$bin_size)

  dat = as.data.frame(cbind(dat$start,apply(dat[,-1],2,function(x){
    ksmooth(dat$start, x, kernel = "normal",
            bandwidth = kernelLen*samplingTime, n.points = length(dat$start))$y
  })))

  # add back all columns
  dat <- reshape2::melt(dat, id.vars = 'V1')
  colnames(dat) <- c('start','sample','value')
  met <- unique(dat_tmp[,-na.omit(match(c('seqnames','start','end','value' ),colnames(dat_tmp) ))])
  dat <- merge(dat,met, by = 'sample')
  dat$end <- dat$start + samplingTime

  divby = 1
  if(mean(seq(start,end, length.out = 5) > 1) > 0.5){
    divby = 1e6
  }

  dat$sample <- factor(dat$sample, levels = colnames(sce))

  if(!is.null(GR)){
    tmp_dat <- as.data.frame(subsetByOverlaps(GR, GRanges(chrom,
                                                          IRanges(start = start,
                                                                  end = end))))
    tmp_dat <- reshape2::melt(tmp_dat[,-c(4:5)], id.vars = c('seqnames','start','end'))
    tmp_dat <- tmp_dat[,c('variable', 'seqnames', 'start', 'end', 'value')]
    colnames(tmp_dat)[1] <- 'sample'
    if(!is.null(colour_on)){
      tmp_dat[[colour_on]] <- tmp_dat$sample
    }
    dat <- rbind(dat[,unique(c('sample', 'seqnames', 'start', 'end', 'value',colour_on))],
                 tmp_dat[,unique(c('sample', 'seqnames', 'start', 'end', 'value',colour_on))])
    dat$sample <- factor(dat$sample, levels = c(colnames(sce), colnames(mcols(GR))))
  }
  if(is.null(ylim)){
    if(ylim_per_sample){

      ds <- split(dat,dat$sample)
      ds <- data.table::rbindlist(lapply(ds, function(x){
        ylim = round(quantile(x$value[is.finite(x$value)],c(.005,.99)),digits = 1)
        x$value[x$value > max(ylim)] <- max(ylim)
        x$value[x$value < min(ylim)] <- min(ylim)
        ylim = NULL
        x
      }))

    } else {
      ylim = quantile(dat$value[is.finite(dat$value)],c(.005,.995))
      dat$value[dat$value > max(ylim)] <- max(ylim)
      dat$value[dat$value < min(ylim)] <- min(ylim)
    }

  }

  if(!is.null(annot_GRlist)){
    dat$value[dat$value < 0] <- 0
  }
  dat$sample <- factor(dat$sample, levels = colnames(sce))
  lbdf <- data.table::rbindlist(lapply(split(dat,dat$sample),function(x){
    data.frame('sample' = unique(x$sample), 'start' = start+((end-start)*0.01), 'value' = floor_dec(max(x$value)))}))
  lbdf$sample <- factor(lbdf$sample, levels = colnames(sce))

  p <- ggplot2::ggplot(dat, ggplot2::aes_string(x = 'start', y = 'value', fill = colour_on)) +
    ggplot2::geom_col(width = median(dat$end - dat$start)) +
    ggplot2::labs(y = legend_label, x= paste0(chrom, ifelse(divby != 1, ' (Mb)',' (bp)'))) +
    ggplot2::theme(aspect.ratio = plot_ratio) +
    ggplot2::geom_hline(yintercept = hline) +
    ggplot2::scale_x_continuous(expand = c(0,0), limits = c(start,end),
                                breaks = seq(start,end, length.out = 5),
                                labels = seq(start,end, length.out = 5)/divby) +
    ggthemes::scale_fill_tableau('Tableau 20')+
    ggplot2::guides(fill = 'none')+
    ggplot2::theme(strip.text.x  = ggplot2::element_blank()) +
    ggplot2::geom_text(data = lbdf,
                       mapping = ggplot2::aes(x = start, label = sample, y = value), vjust = 1, hjust = 0,
                       colour = 'black', size = 6*0.36,family = 'Helvetica')
  if(!is.null(annot_GRlist)){
    annot = as.data.frame(annot)
    annot$sample <- factor(annot$sample, levels = colnames(sce))

    tmp <- lbdf[,-2]
    tmp <- cbind(tmp,t(sapply(tmp$value/3, function(x){quantile(seq(from = 0, to = x,length.out = 10),c(1/3,2/3))})))
    colnames(tmp)[3:4] <- c('yShallow','yDeep')

    annot <- merge(annot,tmp, by = 'sample')

    if(is.null(annot_col)){
      p <- p + ggplot2::geom_rect(data= annot,
                                  #x = start,
                                  #y = min(ylim),
                                  mapping = aes(xmin = start,
                                                xmax = end,
                                                ymin = -1 * yShallow,
                                                ymax = -1 * yDeep ))
    } else {
      p <- p + ggplot2::geom_rect(data= annot,
                                  #x = start,
                                  fill = annot_col,
                                  #y = min(ylim),
                                  mapping = aes(xmin = start,
                                                xmax = end,
                                                ymin = -1 * yShallow,
                                                ymax = -1 * yDeep ))
    }
  }

  if(!ylim_per_sample){
    p <- p + ggplot2::coord_cartesian(ylim = ylim) + ggplot2::facet_grid(forcats::fct_relevel(sample,colnames(sce)) ~ .)
  } else {
    p <- p + ggplot2::facet_grid(forcats::fct_relevel(sample,colnames(sce)) ~ ., scales = 'free_y')
  }

  p  + ggplot2::scale_y_continuous(breaks = function(x){floor_dec(max(x),level = 2)}) + ggplot2::theme(strip.text = ggplot2::element_blank())


}


#' Plot values of a SCE as a raster of pairwise cosine-similarities
#'
#' `mabid_cosplot()` generates a heatmap
#'
#' @param sce An MAb-ID SCE object
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param zlim Set range of z-axis (= fill intensity)
#' @param clustering_method The agglomeration method to be used. This should be
#'     (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'     "complete", "average", "mcquitty", "median" or "centroid".
#' @param rownames Display sample-names?
#' @param size Set size of the heatmap
#' @param font_size Set overall font size
#' @param annotation_columns a vector of column-names of colData(sce)
#' @param annotation_col a vector of hex-colours to use
#' @return A Ggplot2-object
#' @examples
#' \dontrun{mabid_cosplot(sce,assay_to_use = 'counts')}
#' @export
mabid_cosplot <- function(sce, assay_to_use = 'logcounts', zlim = NULL, clustering_method = "ward.D2",
                          rownames = F,
                          size = NULL, font_size = 10, annotation_columns = NULL, annotation_col = NULL){

  if(!is.null(size)){

    if(!grid::is.unit(size)){
      size <- grid::unit(size,'cm')
    }

  } else {
    size = grid::unit(5,'cm')
  }

  if(is.null(annotation_col)){
    annotation_col <- c("#000000", '#777777', '#FFFFFF',"#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'#3b3b3b')
  }

  cpp <- NULL
  annotation_df <- NULL
  if(is.null(annotation_columns)){
    annotation_col <- NULL
  } else {
    annotation_df <- as.data.frame(SingleCellExperiment::colData(sce)[,annotation_columns, drop=FALSE])
    cpp = apply(annotation_df, 2, function(x){setNames(annotation_col[1:length(unique(x))],unique(x))})
  }

  if(!is.null(annotation_df)){
    if(ncol(annotation_df) == 1){
      annotation_col <- lapply(1:ncol(cpp),function(i){cpp[,i]})
      names(annotation_col) <- colnames(cpp)
    } else {
      annotation_col <- cpp
    }

    annotation_col = lapply(annotation_col, function(cp){
      if(any(is.na(cp))){
        nm <- names(cp)
        cp <- colorRampPalette(unname(cp[!is.na(cp)]))(length(cp))
        names(cp) <- nm
      }
      cp
    })
  }


  pal_redAndBlue <- c('#00334a','#21516a', '#44708b', '#6691ad', '#88b3d0', '#b3d5ee',
                      '#f5f5f5',
                      '#fbc5b9', '#f79580', '#de6c5b', '#bf4438', '#98231d', '#6e0000')

  mat     <- assay(sce, assay_to_use)

  # cosine similarity
  similarity <- qlcMatrix::cosSparse(mat)
  similarity[is.na(similarity)] <- 0
  diag(similarity) <- NA
  cdist <- as.dist(1 - similarity)
  hc <- hclust(cdist, method = clustering_method)

  if(is.null(zlim)){
    zlim <- range(similarity,na.rm = T)
  }

  hmap <- NULL
  if(is.null(annotation_col)){
    hmap = ComplexHeatmap::Heatmap(as.matrix(similarity), height = size,width = size,cluster_rows = hc,cluster_columns = hc,
                                   col = circlize::colorRamp2(breaks = seq(from = zlim[1], to = zlim[2],
                                                                           length.out = length(pal_redAndBlue)),
                                                              colors = pal_redAndBlue),
                                   row_names_gp = grid::gpar(fontsize = font_size),
                                   show_row_names = rownames,
                                   show_column_names = F, show_column_dend = F,border = 'black',na_col = tail(pal_redAndBlue,1),
                                   heatmap_legend_param = list(title_position = 'leftcenter-rot',
                                                               title_gp = grid::gpar(fontsize = font_size),
                                                               labels_gp = grid::gpar(fontsize = font_size),
                                                               border = T,
                                                               legend_height = size,border = 'black',
                                                               title = expression(Cos~theta)))
  } else {
    hmap = ComplexHeatmap::Heatmap(as.matrix(similarity), height = size,width = size,cluster_rows = hc,cluster_columns = hc,
                                   col = circlize::colorRamp2(breaks = seq(from = zlim[1], to = zlim[2],
                                                                           length.out = length(pal_redAndBlue)),
                                                              colors = pal_redAndBlue),
                                   show_row_names = rownames,
                                   row_names_gp = grid::gpar(fontsize = font_size),
                                   show_column_names = F, show_column_dend = F,border = 'black',na_col = tail(pal_redAndBlue,1),
                                   heatmap_legend_param = list(title_position = 'leftcenter-rot',
                                                               title_gp = grid::gpar(fontsize = font_size),
                                                               labels_gp = grid::gpar(fontsize = font_size),
                                                               legend_height = size,border = 'black',
                                                               border = T,
                                                               title = expression(Cos~theta)),
                                   left_annotation = ComplexHeatmap::rowAnnotation(annotation_name_gp = grid::gpar(fontsize = font_size),
                                                                                   df = annotation_df[colnames(mat),],
                                                                                   annotation_name_side ='bottom' ,
                                                                                   col=annotation_col,
                                                                                   gap = grid::unit(0, "points"),
                                                                                   border = T,
                                                                                   annotation_legend_param = list(
                                                                                     border = T,
                                                                                     title_gp = grid::gpar(fontsize = font_size),
                                                                                     labels_gp = grid::gpar(fontsize = font_size))))
  }

  ComplexHeatmap::draw(hmap,align_heatmap_legend = 'heatmap_top')
}



#' Plot values of a SCE as a raster of pairwise correlations
#'
#' `mabid_corplot()` generates a heatmap
#'
#' @param sce An MAb-ID SCE object
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param zlim Set range of z-axis (= fill intensity)
#' @param method a character string indicating which correlation coefficient
#'     (or covariance) is to be computed. One of "spearman" (default), "kendall",
#'     or "pearson": CAN NOT be abbreviated.
#' @param clustering_method The agglomeration method to be used. This should be
#'     (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'     "complete", "average", "mcquitty", "median" or "centroid".
#' @param size Set size of the heatmap
#' @param size Set size of the heatmap
#' @param font_size Set overall font size
#' @param rownames Display sample-names?
#' @param annotation_columns a vector of column-names of colData(sce)
#' @param annotation_col a vector of hex-colours to use
#' @return A Ggplot2-object
#' @examples
#' \dontrun{
#' mabid_corplot(sce,assay_to_use = 'counts')
#' }
#' @export
mabid_corplot <- function(sce, assay_to_use = 'logcounts', zlim = NULL, method = "spearman",clustering_method = 'ward.D2',
                          rownames = F,
                          size = NULL, font_size = 10, annotation_columns = NULL, annotation_col = NULL){

  if(!is.null(size)){

    if(!grid::is.unit(size)){
      size <- grid::unit(size,'cm')
    }

  } else {
    size = grid::unit(5,'cm')
  }

  if(is.null(annotation_col)){
    annotation_col <- c("#000000", '#777777', '#FFFFFF',"#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'#3b3b3b')
  }

  legend_label <- switch(method,
                         'pearson' = expression(italic(r)),
                         'kendall' = expression(tau),
                         'spearman' = expression(rho))

  cpp <- NULL
  annotation_df <- NULL
  if(is.null(annotation_columns)){
    annotation_col <- NULL
  } else {
    annotation_df <- as.data.frame(SingleCellExperiment::colData(sce)[,annotation_columns, drop=FALSE])
    cpp = apply(annotation_df, 2, function(x){setNames(annotation_col[1:length(unique(x))],unique(x))})
  }


  if(!is.null(annotation_df)){
    if(ncol(annotation_df) == 1){
      annotation_col <- lapply(1:ncol(cpp),function(i){cpp[,i]})
      names(annotation_col) <- colnames(cpp)
    } else {
      annotation_col <- cpp
    }

    annotation_col = lapply(annotation_col, function(cp){
      if(any(is.na(cp))){
        nm <- names(cp)
        cp <- colorRampPalette(unname(cp[!is.na(cp)]))(length(cp))
        names(cp) <- nm
      }
      cp
    })
  }

  pal_redAndBlue <- c('#00334a','#21516a', '#44708b', '#6691ad', '#88b3d0', '#b3d5ee',
                      '#f5f5f5',
                      '#fbc5b9', '#f79580', '#de6c5b', '#bf4438', '#98231d', '#6e0000')

  mat     <- as.matrix(assay(sce, assay_to_use))

  # cosine similarity
  correlation <- cor(mat, method = method, use = "all.obs")
  correlation[is.na(correlation)] <- 0
  diag(correlation) <- NA
  cdist <- as.dist(1 - correlation)
  hc <- hclust(cdist, method = clustering_method)

  if(is.null(zlim)){
    zlim <- range(correlation,na.rm = T)
  }

  hmap <- NULL
  if(is.null(annotation_col)){
    hmap = ComplexHeatmap::Heatmap(as.matrix(correlation), height = size,width = size,cluster_rows = hc,cluster_columns = hc,
                                   col = circlize::colorRamp2(breaks = seq(from = zlim[1], to = zlim[2],
                                                                           length.out = length(pal_redAndBlue)),
                                                              colors = pal_redAndBlue),
                                   row_names_gp = grid::gpar(fontsize = font_size),
                                   show_row_names = rownames,
                                   show_column_names = F, show_column_dend = F,border = 'black',na_col = tail(pal_redAndBlue,1),
                                   heatmap_legend_param = list(title_position = 'leftcenter-rot',
                                                               title_gp = grid::gpar(fontsize = font_size),
                                                               border = T,
                                                               labels_gp = grid::gpar(fontsize = font_size),
                                                               legend_height = size,border = 'black',
                                                               title = legend_label))
  } else {
    hmap = ComplexHeatmap::Heatmap(as.matrix(correlation), height = size,width = size,cluster_rows = hc,cluster_columns = hc,
                                   col = circlize::colorRamp2(breaks = seq(from = zlim[1], to = zlim[2],
                                                                           length.out = length(pal_redAndBlue)),
                                                              colors = pal_redAndBlue),
                                   show_row_names = rownames,
                                   row_names_gp = grid::gpar(fontsize = font_size),
                                   show_column_names = F, show_column_dend = F,border = 'black',na_col = tail(pal_redAndBlue,1),
                                   heatmap_legend_param = list(title_position = 'leftcenter-rot',
                                                               title_gp = grid::gpar(fontsize = font_size),
                                                               labels_gp = grid::gpar(fontsize = font_size),
                                                               border = T,
                                                               legend_height = size,border = 'black',
                                                               title = legend_label),
                                   left_annotation = ComplexHeatmap::rowAnnotation(annotation_name_gp = grid::gpar(fontsize = font_size),
                                                                                   df = annotation_df[colnames(mat),],
                                                                                   border = T,
                                                                                   gap = grid::unit(0, "points"),
                                                                                   annotation_name_side ='bottom',
                                                                                   col=annotation_col,
                                                                                   annotation_legend_param = list(
                                                                                     border = T,
                                                                                     title_gp = grid::gpar(fontsize = font_size),
                                                                                     labels_gp = grid::gpar(fontsize = font_size))))
  }

  ComplexHeatmap::draw(hmap,align_heatmap_legend = 'heatmap_top')
}

#' Plot values of a ROI from a SCE as rasterplot
#'
#' `mabid_rasterplot()` generates a heatmap
#'
#' Pretty pretty!
#'
#' @param sce An MAb-ID SCE object
#' @param chrom Chromosome of interest
#' @param start Start of region
#' @param end End of region
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param clustering_method The agglomeration method to be used. This should be
#'     one of "ward.D", "ward.D2", "single",
#'     "complete", "average", "mcquitty", "median" or "centroid".
#' @param clustering_distance the distance measure to be used. This must be one of
#'     "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param divpal Use a divergent colourmap?
#' @param zlim Set range of z-axis (= fill intensity)
#' @param asp The aspect-ratio of the plot
#' @param size Set size of the heatmap
#' @param font_size Set overall font size
#' @param rownames Display sample-names?
#' @param annotation_columns a vector of column-names of colData(sce)
#' @param annotation_col a vector of hex-colours to use
#' @return A Ggplot2-object
#' @examples
#' \dontrun{
#' mabid_rasterplot(sce, end = 249250621, assay_to_use = 'counts')
#' }
#' @export
#'
#'
mabid_rasterplot <- function(sce, chrom = 'chr1', start = 0, end = 5000000,
                             assay_to_use = 'logcounts',
                             clustering_method = 'ward.D2',
                             clustering_distance = "canberra",
                             divpal = F, zlim = NULL, asp = 1.618034, size = NULL, font_size = 10,
                             rownames = F,
                             annotation_columns = NULL, annotation_col = NULL){

  legend_label <- switch(assay_to_use,
                         'OE_control_log2' = expression(Log[2]~(frac(counts, control))),
                         'OE_mappability_log2' = expression(Log[2]~(frac(counts, mappability))),
                         'counts' = 'counts','normcounts' = 'normalised signal',
                         'logcounts' = expression(ln~(CPM)),
                         assay_to_use)

  ############################################################################## set colours
  pal_redAndBlue <- c('#00334a', '#21516a', '#44708b', '#6691ad', '#88b3d0', '#b3d5ee',
                      '#f5f5f5',
                      '#fbc5b9', '#f79580', '#de6c5b', '#bf4438', '#98231d', '#6e0000')
  pal_blue <- pal_redAndBlue[7:1]
  pal_grey <- colorspace::desaturate(pal_blue)

  CP = pal_grey#inlmisc::GetColors(50,'iridescent')
  if(divpal){
    CP = pal_redAndBlue# inlmisc::GetColors(50,'PRGn')
  }

  ############################################################################## get annotations
  if(is.null(annotation_col)){
    annotation_col <- c("#000000", '#777777', '#FFFFFF',"#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'#3b3b3b')
  }
  cpp <- NULL
  annotation_df <- NULL
  if(is.null(annotation_columns)){
    annotation_col <- NULL
  } else {
    annotation_df <- as.data.frame(SingleCellExperiment::colData(sce)[,annotation_columns, drop=FALSE])
    cpp = apply(annotation_df, 2, function(x){setNames(annotation_col[1:length(unique(x))],unique(x))})
  }

  if(!is.null(annotation_df)){
    if(ncol(annotation_df) == 1){
      annotation_col <- lapply(1:ncol(cpp),function(i){cpp[,i]})
      names(annotation_col) <- colnames(cpp)
    } else {
      annotation_col <- cpp
    }

    annotation_col = lapply(annotation_col, function(cp){
      if(any(is.na(cp))){
        nm <- names(cp)
        cp <- colorRampPalette(unname(cp[!is.na(cp)]))(length(cp))
        names(cp) <- nm
      }
      unlist(cp)
    })
  }



  ############################################################################## get data
  dat <- mabid_get_track(sce, chrom, start, end, assay_to_use,metadata = F)
  mat <- tidyr::spread(dat[,c(2,4,5)], key = 'sample',value = 'value')
  mat <- mat[order(mat[,1]),]
  mat <- as.matrix(mat); rownames(mat) <- mat[,1]; mat <- mat[,-1];mat <- t(mat)

  if(is.null(zlim)){
    zlim = range(dat$value)
  }
  dat$value[dat$value > max(zlim)] <- max(zlim)
  dat$value[dat$value < min(zlim)] <- min(zlim)

  ############################################################################## get position-markers
  all_column_breaks <- as.numeric(colnames(mat))
  breaks <- c(start,end)

  divby = 1
  if((breaks[2] - breaks[1]) > 1e6){
    divby = 1e6
  }
  breaks <- round(breaks/divby,1)
  labels <- sub("\\.0+$", "", as.character(breaks))
  column_labels <- rep('',length(all_column_breaks))
  column_labels[range(seq_along(column_labels))] <- labels

  ############################################################################## draw heatmap
  if(!is.null(size)){

    if(!grid::is.unit(size)){
      size <- grid::unit(size,'cm')
    }

  } else {
    size = grid::unit(5,'cm')
  }

  hmap <- NULL
  if(is.null(annotation_col)){
    hmap = ComplexHeatmap::Heatmap(mat,height = size,width = size*asp,cluster_columns = F,column_labels = column_labels,
                                   col = circlize::colorRamp2(breaks = seq(from = zlim[1], to = zlim[2],
                                                                           length.out = length(CP)),
                                                              colors = CP),column_title_side = 'bottom',
                                   clustering_distance_rows = clustering_distance,
                                   column_title = paste0(chrom, ifelse(divby != 1, ' (Mb)',' (bp)')),
                                   row_names_gp = grid::gpar(fontsize = font_size),clustering_method_rows = clustering_method,
                                   show_column_dend = F,border = 'black',na_col = tail(pal_redAndBlue,1),column_names_centered = T,
                                   column_names_gp = grid::gpar(fontsize = font_size),column_names_rot = 0,
                                   column_title_gp = grid::gpar(fontsize = font_size),show_row_names = rownames,
                                   heatmap_legend_param = list(title_position = 'leftcenter-rot',
                                                               title_gp = grid::gpar(fontsize = font_size),
                                                               labels_gp = grid::gpar(fontsize = font_size),
                                                               legend_height = size,border = 'black',
                                                               title = legend_label))
  } else {
    hmap = ComplexHeatmap::Heatmap(mat,height = size,width = size*asp,cluster_columns = F,column_labels = column_labels,
                                   col = circlize::colorRamp2(breaks = seq(from = zlim[1], to = zlim[2],
                                                                           length.out = length(CP)),
                                                              colors = CP),column_title_side = 'bottom',
                                   clustering_distance_rows = clustering_distance,
                                   column_title = paste0(chrom, ifelse(divby != 1, ' (Mb)',' (bp)')),
                                   row_names_gp = grid::gpar(fontsize = font_size),clustering_method_rows = clustering_method,
                                   show_column_dend = F,border = 'black',na_col = tail(pal_redAndBlue,1),column_names_centered = T,
                                   column_names_gp = grid::gpar(fontsize = font_size),column_names_rot = 0,
                                   column_title_gp = grid::gpar(fontsize = font_size),show_row_names = rownames,
                                   heatmap_legend_param = list(title_position = 'leftcenter-rot',
                                                               title_gp = grid::gpar(fontsize = font_size),
                                                               labels_gp = grid::gpar(fontsize = font_size),
                                                               legend_height = size,border = 'black',
                                                               title = legend_label),
                                   left_annotation = ComplexHeatmap::rowAnnotation(annotation_name_gp = grid::gpar(fontsize = font_size),
                                                                                   df = annotation_df[rownames(mat),],
                                                                                   annotation_name_side ='bottom',
                                                                                   col=annotation_col,
                                                                                   border = T,
                                                                                   gap = grid::unit(0, "points"),
                                                                                   annotation_legend_param = list(
                                                                                     border = T,
                                                                                     title_gp = grid::gpar(fontsize = font_size),
                                                                                     labels_gp = grid::gpar(fontsize = font_size))
                                   ))
  }
  ComplexHeatmap::draw(hmap,align_heatmap_legend = 'heatmap_top')

}





#' Plot values of a ROI from a SCE as raster, with a population-average on top
#'
#' `mabid_populationplot()` generates a heatmap
#'
#' Pretty pretty!
#'
#' @param sce An MAb-ID SCE object
#' @param chrom Chromosome of interest
#' @param start Start of region
#' @param end End of region
#' @param group_by Select one (!!!) colData-column to split on
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param ylim Set range of y-axis
#' @param zlim Set range of z-axis (= fill intensity)
#' @param hline Add horizontal line to upper plot
#' @param vline Add vertical line to both plots
#' @param line_colour Set line colour
#' @param line_type Set line type
#' @param plot_ratio Set height-ratios of upper versus lower plot (default is 1/10)
#' @param FUN Function to aggregate over. Set to mean
#' @param title Title
#' @param subtitle Subtitle
#' @param divpal Use a divergent colourmap?
#' @return A Ggplot2-object
#' @examples
#' \dontrun{
#' mabid_populationplot(sce, end = 249250621, assay_to_use = 'counts')
#' }
#' @export
#'
#'
mabid_populationplot <- function(sce,
                             chrom = 'chr1',
                             start = 20e6,
                             end = 30e6,
                             group_by = NULL,# max 1!
                             assay_to_use = 'logcounts',
                             ylim = NULL,
                             zlim = NULL,
                             hline = 0,
                             vline = NULL,
                             line_colour = 'black',
                             line_type = 1,
                             plot_ratio = 1/10,
                             FUN = mean,
                             title = NULL,
                             subtitle = NULL,
                             divpal = F){
  require(patchwork)
  legend_label <- switch(assay_to_use,
                         'OE_control_log2' = expression(Log[2]~(frac(counts, control))),
                         'OE_mappability_log2' = expression(Log[2]~(frac(counts, mappability))),
                         'counts' = 'counts',
                         'normcounts' = 'normalised signal',
                         'logcounts' = expression(ln~(CPM)),
                         assay_to_use)
  ############################################################################## set colours
  pal_redAndBlue <- c('#00334a', '#21516a', '#44708b', '#6691ad', '#88b3d0', '#b3d5ee',
                      '#f5f5f5',
                      '#fbc5b9', '#f79580', '#de6c5b', '#bf4438', '#98231d', '#6e0000')
  pal_blue <- pal_redAndBlue[7:1]
  pal_grey <- c("#F5F5F5", "#D0D0D0", "#ADADAD", "#8C8C8C", "#6B6B6B", "#4C4C4C",
                "#2F2F2F") #colorspace::desaturate(pal_blue)

  CP = pal_grey#inlmisc::GetColors(50,'iridescent')
  if(divpal){
    CP = pal_redAndBlue# inlmisc::GetColors(50,'PRGn')
  }

  ############################################################################## get data
  colData(sce) <- colData(sce)[,c('bin_size','depth',group_by)]

  dat <- data.table::data.table(mabidR::mabid_get_track(sce,chrom,start,
                                                        end+(sce$bin_size[1]*2),assay_to_use,metadata = T))
  dat <- dat[ , FUN(value),
              by = eval(paste0(c(c('sample','start','end','depth'),group_by),collapse = ','))]

  if(is.null(group_by)){
    dat$group <- 'no_group'
    dat <- dat[,c(1,2,3,4,6,5)]
    colnames(dat)[6] <- c('value')
  } else {
    colnames(dat)[5:6] <- c('group','value')
  }
  dat_population <- dat[,FUN(value), by = eval(paste0(c('start','end','group'),collapse = ','))]
  colnames(dat_population)[4] <- c('value')

  divby = 1
  if(mean(seq(start,end, length.out = 5) > 1) > 0.5){
    divby = 1e6
  }

  if(is.null(ylim)){
    ylim = quantile(dat_population$value,c(.005,.995))
  }
  if(is.null(zlim)){
    zlim = quantile(dat$value,c(.005,.995))
  }

  dat_population$value[dat_population$value > max(ylim)] <- max(ylim)
  dat_population$value[dat_population$value < min(ylim)] <- min(ylim)
  dat$value[dat$value > max(zlim)] <- max(zlim)
  dat$value[dat$value < min(zlim)] <- min(zlim)

  p1 <- ggplot2::ggplot(dat_population, ggplot2::aes_string(x = 'start', y = 'value', col = 'group')) +
    ggplot2::geom_line() +
    ggplot2::labs(y = NULL, x= NULL, colour = NULL, title = title, subtitle = subtitle) +
    ggplot2::geom_hline(yintercept = hline, colour = line_colour,lty = line_type) +
    ggplot2::geom_vline(xintercept = vline, colour = line_colour,lty = line_type) +
    ggplot2::scale_x_continuous(expand = c(0,0), breaks = seq(start,end, length.out = 5),
                                labels = seq(start,end, length.out = 5)/divby) +
    # ggplot2::scale_color_manual(values = colorRampPalette(erithacus::colour_blind_palette())(length(unique(dat_population$group))))+
    ggplot2::coord_cartesian(ylim = ylim, xlim = c(start,end)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = 'black'),
      axis.line = ggplot2::element_blank(),

      legend.key = ggplot2::element_rect(fill = NA),

      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = 'black', size = 1),
      panel.grid = ggplot2::element_blank(),

      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = 'black'),

      text = ggplot2::element_text(color = "black")
    )

  if(is.null(group_by)){
    p1 <- p1 + ggplot2::guides(colour = 'none')
  }

  # vergeet niet te sorteren!


  dat <- data.table::rbindlist(lapply(split(dat,dat$group), function(x){
    srt <- unique(x[,c('sample','depth')])
    srt$order <- factor(order(srt$depth,decreasing = T))
    x <- merge(x,srt[,-2])
    x
  }))


  p2 <- ggplot2::ggplot(dat,  ggplot2::aes_string(x = 'start', y = 'order', fill = 'value')) +
    ggrastr::geom_tile_rast()+
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0), breaks = seq(start,end, length.out = 5),
                                labels = seq(start,end, length.out = 5)/divby) +
    ggplot2::scale_fill_gradientn(colours = CP,)+
    ggplot2::geom_vline(xintercept = vline, colour = line_colour,lty = line_type) +
    ggplot2::labs(fill = legend_label,y = NULL,x= paste0(chrom, ifelse(divby != 1, ' (Mb)',' (bp)'))) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = 'black'),
      axis.line = ggplot2::element_blank(),

      legend.key = ggplot2::element_rect(fill = NA),

      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = 'black', size = 1),
      panel.grid = ggplot2::element_blank(),

      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = 'black'),

      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),

      text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::coord_cartesian(xlim = c(start,end)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(frame.linewidth = 1,frame.colour = 'black',ticks.colour = 'black'))

  if(!is.null(group_by)){
    p2 <- p2 + ggplot2::facet_grid(rows = 'group',space = 'free_y', scales = 'free_y')
  }

  p1 / p2 + patchwork::plot_layout(ncol = 1,widths = 1,guides = "collect",
                                   heights = c(10*plot_ratio,
                                                        10*(1-plot_ratio)))
}


#' Performd SVD and plot results
#'
#' `mabid_SVD()` returns a ggplot2-object
#'
#' Pretty pretty!
#'
#' @param sce An MAb-ID SCE object
#' @param svx Which SV on the x-axis?
#' @param svy Which SV on the y-axis?
#' @param assay_to_use Use data from one of the assays the sce (default: `logcounts`). Look in `names(assays(sce))`
#' @param colour_by Use a column of the output of `mabid_get_track(sce)` (default: 'sample')
#' @param shape_by Use a column of the output of `mabid_get_track(sce)` (default: 'sample')
#' @return A Ggplot2-object
#' @examples
#' \dontrun{
#' mabid_SVD(sce_SBC001)
#' }
#' @export
mabid_SVD <- function(sce, svx = 'SV_1', svy = 'SV_2', colour_by = NULL, shape_by = NULL,
                      assay_to_use = 'logcounts'){

  x <- assay(sce, assay_to_use)
  x_ <- svd(as.matrix(x))
  x_ <- x_$v
  rownames(x_) <- colnames(x)
  x_ <- as.data.frame(x_)
  colnames(x_) <- gsub('V', 'SV_', colnames(x_))
  x_$id <- rownames(x_)

  if (!is.null(colour_by)) {
    info <- setNames(colData(sce)[, colour_by], colnames(sce))
    x_$kleur <- unname(info[x_$id])
  } else {
    x_$kleur <- x_[, 1]
  }

  if (!is.null(shape_by)) {
    info <- setNames(colData(sce)[, shape_by], colnames(sce))
    x_$vorm <- unname(info[x_$id])
  } else {
    x_$vorm <- x_[, 1]
  }
  if (length(unique(x_$vorm)) > 5) {
    stop('shape_by must have <= 5 unique values')
  }
  scale_vorm <- setNames(c(16, 15, 17, 18, 19), unique(x_$vorm))
  scale_vorm <- scale_vorm[!is.na(names(scale_vorm))]

  x_$svx <- x_[, svx]
  x_$svy <- x_[, svy]

  ggplot2::ggplot(x_, ggplot2::aes(
    x = svx,
    y = svy,
    col = kleur,
    shape = vorm
  )) +
    # erithacus::erithacus_theme(ratio = 'equal') +
    ggplot2::labs(x = svx, y = svy) +
    #ggplot2::scale_colour_manual(values = erithacus::colour_blind_palette()) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_shape_manual(values = scale_vorm) +
    ggplot2::labs(
      x = gsub('_', ' ', svx),
      y = gsub('_', ' ', svy),
      col = colour_by,
      shape = shape_by
    )
}







round.off <- function (x, digits=0) {
  posneg = sign(x)
  z = trunc(abs(x) * 10 ^ (digits + 1)) / 10
  z = floor(z * posneg + 0.5) / 10 ^ digits
  return(z)
}

#' Plot values of a ROI from a SCE as mirrors
#'
#' `mabid_mirrorplot()` returns generates mirrorplots.
#'
#' Pretty pretty!
#'
#' @param sce An MAb-ID SCE object
#' @param chrom Chromosome of interest
#' @param start Start of region
#' @param end End of region
#' @param mirror_list A named list of named two-length vectors of sample names.
#' @param tip Show a tip on how to save the ggplot.
#' @param panel_width How wide should the
#' @param smooth_kernel Average over N bins
#' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' @param GR add bigwig-like tracks.
#' @param ylim Set range of y-axis
#' @param plot_ratio Set aspect ratio (default: 1/1.618)
#' @return A Ggplot2-object
#' @examples
#' \dontrun{
#'
#' mirror_list <- list("H3K36me3" = c('single' = "H3K36me3: singlet",'combined' = "H3K36me3: combined"),
#' 'H3K27me3' = c('single' = "H3K27me3: singlet", 'combined' = "H3K27me3: combined"),
#' 'PolII S5P' = c('single' = 'RNA PolII CTD S5P: singlet','combined' = 'RNA PolII CTD S5P: combined'),
#' 'H3' = c('single' = "H3: singlet", 'combined' = "H3: combined"))
#' mabid_mirrorplot(sce_SBC001,mirror_list = mirror_list,
#'                            chrom = 'chr1', start = 0, end = 100e6,
#'                            assay_to_use = 'logcounts')
#' }
#' @export
mabid_mirrorplot <- function(sce, chrom = 'chr1', start = 0, end = 50000000, mirror_list, tip = T, panel_width = 1,smooth_kernel = 1,
                             assay_to_use = 'logcounts', GR = NULL, ylim = NULL, plot_ratio = 1/1.618){

  if(!length(unique(names(mirror_list))) == length(names(mirror_list))){
    stop('Names in mirror_list should be unique!')
  }

  if(!all(sapply(mirror_list, length) == 2)){
    stop('Not all entries in mirror_list have two sample-names!')
  }

  legend_label <- switch(assay_to_use,
                         'OE_control_log2' = expression(Log[2]~(frac(counts, control))),
                         'OE_mappability_log2' = expression(Log[2]~(frac(counts, mappability))),
                         'counts' = 'counts','normcounts' = 'normalised signal',
                         'logcounts' = expression(ln~(CPM)),
                         assay_to_use)

  dat <- mabid_get_track(sce, chrom, start, end, assay_to_use)
  divby = 1
  if(mean(seq(start,end, length.out = 5) > 1) > 0.5){
    divby = 1e6
  }


  dat$sample <- factor(dat$sample, levels = colnames(sce))

  if(!is.null(GR)){
    tmp_dat <- as.data.frame(subsetByOverlaps(GR, GRanges(chrom,
                                                          IRanges(start = start,
                                                                  end = end))))
    tmp_dat <- reshape2::melt(tmp_dat[,-c(4:5)], id.vars = c('seqnames','start','end'))
    tmp_dat <- tmp_dat[,c('variable', 'seqnames', 'start', 'end', 'value')]
    colnames(tmp_dat)[1] <- 'sample'
    dat <- rbind(dat[,unique(c('sample', 'seqnames', 'start', 'end', 'value'))],
                 tmp_dat[,unique(c('sample', 'seqnames', 'start', 'end', 'value'))])
    dat$sample <- factor(dat$sample, levels = c(colnames(sce), colnames(mcols(GR))))
  }

  if(!all(unname(unlist(lapply(mirror_list, unique))) %in% dat$sample)){
    stop('Not all sample-names in mirror_list are found.')
  }


  # do the magic mirror-thing
  track <- data.table::rbindlist(lapply(seq_along(mirror_list), function(i){
    track_top <- dat[dat$sample == mirror_list[[i]][1],c("start","end","value")]
    track_bot <- dat[dat$sample == mirror_list[[i]][2],c("start","end","value")]

    track_top$value[track_top$value < 0] <- 0
    track_bot$value[track_bot$value < 0] <- 0

    track_bot$value <- track_bot$value * -1

    track_top$sample <- names(mirror_list[[i]][1])
    track_bot$sample <- names(mirror_list[[i]][2])

    track_top$placing <- 'top'
    track_bot$placing <- 'bottom'

    track <- rbind(track_top,track_bot)
    track$mirr <- names(mirror_list)[i]
    track
  }))


  if(is.null(ylim)){
    ylim = round.off(quantile(abs(track$value),.995),digits = 1)
  }
  if(length(ylim) == 1) {
    ylim = c(ylim,ylim)
  }
  ylim <- unname(ylim)
  track$value[track$value > ylim[1]] <- ylim[1]
  track$value[track$value < (-1*ylim[2])] <- (-1*ylim[2])

  ann <- unique(track[,c('sample','mirr','placing')])
  ann$start = start; ann$end = start
  ann$value <- ifelse(ann$placing == 'top',ylim[2], -1*ylim[1])
  track$mirr <- factor(track$mirr, levels = names(mirror_list))
  ggplot2::ggplot(track, ggplot2::aes(x = start, y = value, fill = sample)) +
    ggplot2::facet_grid(mirr ~ .) +
    ggseas::stat_rollapplyr(geom = 'area',width = smooth_kernel,align = 'center',FUN = mean) +
    ggplot2::labs(y = legend_label, x= paste0(chrom, ifelse(divby != 1, ' (Mb)',' (bp)'))) +
    # erithacus::theme_npg(box = T,panel_width = panel_width,asp = plot_ratio,tip = tip)+
    ggplot2::guides(fill = "none") +
    ggplot2::geom_hline(yintercept = 0,colour = 'black', size = 0.25) +
    ggplot2::scale_x_continuous(expand = c(0,0), limits = c(start,end),
                                breaks = seq(start,end, length.out = 5),
                                labels = seq(start,end, length.out = 5)/divby) +
    ggplot2::coord_cartesian(ylim = c(-1*ylim[1],ylim[2])) +
    ggplot2::scale_y_continuous(expand = c(0,0),
                                breaks = c(-1 * ylim[1],0,ylim[2]),
                                labels = c(ylim[1],0,ylim[2])) +

    ggplot2::geom_text(ann, mapping = ggplot2::aes(label = sample,vjust = ifelse(value > 0, 1.5,-.5)),
              size = grid::unit(6,'pt')*0.36,hjust = -0.1) +
    ggplot2::scale_fill_manual(values = rep(c('#737373','#525252'),sum(ann$placing == 'top')))
}
