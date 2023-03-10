% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_mabid.R
\name{mabid_cosplot}
\alias{mabid_cosplot}
\title{Plot values of a SCE as a raster of pairwise cosine-similarities}
\usage{
mabid_cosplot(
  sce,
  assay_to_use = "logcounts",
  zlim = NULL,
  clustering_method = "ward.D2",
  rownames = F,
  size = NULL,
  font_size = 10,
  annotation_columns = NULL,
  annotation_col = NULL
)
}
\arguments{
\item{sce}{An MAb-ID SCE object}

\item{assay_to_use}{Use data from one of the assays the sce. Look in `names(assays(sce))`}

\item{zlim}{Set range of z-axis (= fill intensity)}

\item{clustering_method}{The agglomeration method to be used. This should be
(an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
"complete", "average", "mcquitty", "median" or "centroid".}

\item{rownames}{Display sample-names?}

\item{size}{Set size of the heatmap}

\item{font_size}{Set overall font size}

\item{annotation_columns}{a vector of column-names of colData(sce)}

\item{annotation_col}{a vector of hex-colours to use}
}
\value{
A Ggplot2-object
}
\description{
`mabid_cosplot()` generates a heatmap
}
\examples{
\dontrun{mabid_cosplot(sce,assay_to_use = 'counts')}
}
