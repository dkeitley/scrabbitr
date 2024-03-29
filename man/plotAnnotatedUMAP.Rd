% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotAnnotatedUMAP}
\alias{plotAnnotatedUMAP}
\title{Plot a UMAP with annotated labels
Used to show predicted cell type annotations, where the number of factors is
large.}
\usage{
plotAnnotatedUMAP(
  sce,
  colour_by = "predicted_celltype",
  ncell_filt = 5,
  size = 1,
  label_force = 15,
  text_size = 2,
  line_size = 0.5,
  palette = NULL
)
}
\description{
Plot a UMAP with annotated labels
Used to show predicted cell type annotations, where the number of factors is
large.
}
