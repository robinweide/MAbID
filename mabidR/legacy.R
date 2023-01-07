#' #' Plot values of a ROI from a SCE as rasterplot
#' #'
#' #' `mabid_rasterplot()` generates a heatmap
#' #'
#' #' Pretty pretty!
#' #'
#' #' @param sce An MAb-ID SCE object
#' #' @param chrom Chromosome of interest
#' #' @param start Start of region
#' #' @param end End of region
#' #' @param assay_to_use Use data from one of the assays the sce. Look in `names(assays(sce))`
#' #' @param na_col Set NA-values to a specific colour (default: black)
#' #' @param zlim Set range of z-axis (= fill intensity)
#' #' @param hline Add a horizontal line at a y-axis value
#' #' @param plot_ratio Set aspect ratio (default: 1/10)
#' #' @return A Ggplot2-object
#' #' @examples
#' #' \dontrun{
#' #' mabid_rasterplot(sce, end = 249250621, assay_to_use = 'counts')
#' #' }
#' #' @export
#' mabid_rasterplot <- function(sce, chrom = 'chr1', start = 0, end = 5000000, raster = T, assay_to_use = 'OE_control_log2', divpal = F,na_col = NULL, zlim = NULL, hline = 0,plot_ratio = 1/10,round_label = 2){
#'
#'   legend_label <- switch(assay_to_use,
#'                          'OE_control_log2' = expression(Log[2]~(frac(counts, control))),
#'                          'OE_mappability_log2' = expression(Log[2]~(frac(counts, mappability))),
#'                          'counts' = 'counts',
#'                          'logcounts' = expression(ln~(CPM)),
#'                          assay_to_use)
#'
#'   pal_redAndBlue <- c('#00334a', '#21516a', '#44708b', '#6691ad', '#88b3d0', '#b3d5ee',
#'                       '#f5f5f5',
#'                       '#fbc5b9', '#f79580', '#de6c5b', '#bf4438', '#98231d', '#6e0000')
#'   pal_blue <- pal_redAndBlue[7:1]
#'   pal_grey <- colorspace::desaturate(pal_blue)
#'
#'   CP = pal_blue#inlmisc::GetColors(50,'iridescent')
#'   if(divpal){
#'     CP = pal_redAndBlue# inlmisc::GetColors(50,'PRGn')
#'   }
#'
#'   if(is.null(na_col)){
#'     na_col <- 'white' #attr(CP,"nan")
#'     if(is.na(na_col)){
#'       na_col <- 'white'
#'     }
#'   }
#'   dat <- mabid_get_track(sce, chrom, start, end, assay_to_use)
#'   divby = 1
#'   if(mean(seq(start,end, length.out = 5) > 1) > 0.5){
#'     divby = 1e6
#'   }
#'   if(is.null(zlim)){
#'     zlim = range(dat$value)
#'   }
#'   dat$mid <- apply(dat[,c('start','end')], 1, mean)
#'   dat$value[dat$value > max(zlim)] <- max(zlim)
#'   dat$value[dat$value < min(zlim)] <- min(zlim)
#'   dat$sample <- factor(dat$sample, levels = rev(colnames(sce)))
#'   bbb <- round(seq(zlim[1],zlim[2],length.out = 3),round_label)
#'   p = ggplot2::ggplot(dat, ggplot2::aes(x = start, y = sample, fill = value)) +
#'     ggplot2::labs( y = NULL, fill = assay_to_use, x= paste0(chrom, ifelse(divby != 1, ' (Mb)',' (bp)'))) +
#'     erithacus::erithacus_theme(landscape = T)+
#'     #guides(fill = F) +
#'     ggplot2::theme(aspect.ratio = plot_ratio, axis.ticks.y = ggplot2::element_blank()) +
#'     ggplot2::geom_hline(yintercept = hline) +
#'     ggplot2::scale_x_continuous(expand = c(0,0), limits = c(start,end), breaks = seq(start,end, length.out = 5), labels = seq(start,end, length.out = 5)/divby) +
#'     ggplot2::scale_fill_gradientn(breaks = bbb,
#'                                   colours = CP,limits = zlim, na.value = na_col,
#'                                   guide = ggplot2::guide_colourbar(title = legend_label, title.hjust = 0.5,frame.linewidth = 1,
#'                                                                    title.position = 'top',frame.colour = 'black',nbin = 100,
#'                                                                    direction = "horizontal",ticks = F,draw.ulim = T,draw.llim = T)) +
#'     ggplot2::scale_y_discrete(expand = c(0,0)) +
#'     ggplot2::theme(legend.position = 'right', legend.justification = 'top')
#'   if(raster){
#'     p = p + ggrastr::geom_tile_rast()
#'   } else {
#'     p = p + ggplot2::geom_tile()
#'   }
#'   p
#' }
