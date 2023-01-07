#' clean_concentration_annotations
#'
#' `clean_concentration_annotations()`
#'
#' Standarizes concentration-values
#'
#' @param x A character-vector with a hodgepodge of numerical annotations (e.g., "10", "1_in_100", "0,5"))
#' @return A vector of factors
#' @export
clean_concentration_annotations <- function(x){
  x <- gsub(',','.',x)
  y = reshape2::colsplit(x,pattern = "_in_", letters)[,1:2]
  if(length(x) == 1){
    y = rbind(y,y)
  }
  y = apply(y, 2, as.numeric)
  y[,2][is.na(y[,2])] <- 1
  y = y[,1]/y[,2]
  if(length(x) == 1){
    y = y[1]
  }
  factor(y, levels = sort(unique(y)))
}

#' mabid_annotate
#'
#' `mabid_annotate()` lets you add a named vector to colData
#'
#' Easiest why to add some extra annotations
#'
#' @param sce The mabid-sce
#' @param anno_vec A named vector
#' @param name The name of the new annotation-colun
#' @param where Can be set to match to a column-name in coldata, instead of colnames(sce)
#' @return a new sce
#' @export
mabid_annotate <- function(sce, anno_vec, name, where = NULL){
  y = NULL
  if(is.null(where)){
    y <- match(colnames(sce),names(anno_vec))
  } else {
    if(!where %in% colnames(colData(sce))){
      stop('where-argument not a column-name in colData!')
    }
    y <- match(sce[[where]],names(anno_vec))
  }
  colData(sce)[[name]] <- anno_vec[y]
  sce
}

#' mabid_load_sample_sheet
#'
#' Load and parse a sample-sheet
#'
#' @param path The location of the sample-sheet
#' @return a data.table
#' @export
mabid_load_sample_sheet <- function(path) {
  SL <- data.table::fread(data.table = F,path)
  SL$illumina_index = paste0('index_',
                             stringi::stri_pad(SL$illumina_index, width = 2, pad = '0'))
  SL$ABBC_barcode = paste0('ABBC_',
                           stringi::stri_pad(SL$ABBC_barcode, width = 3, pad = '0'))
  SL$SBC_barcode = paste0('SBC_', stringi::stri_pad(SL$SBC_barcode, width = 3, pad = '0'))
  rownames(SL) <- apply(SL[, c("illumina_index","SBC_barcode","ABBC_barcode")], 1, paste0, collapse = '.')
  SL$ID <- rownames(SL)
  SL
}
