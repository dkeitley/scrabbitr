
#' Average expression values within each neighbourhood
#' @importFrom Matrix Diagonal
#' @param X Matrix of expression values
#' @param mask Indicator matrix specifying cell groupings
calcNhoodMean <- function(X, mask) {
  ncells <- Matrix::colSums(mask)
  ncells_inv <- 1/ncells
  ncells_mat <- Matrix::Diagonal(x=ncells_inv)

  nhood_sum <- X %*% mask
  nhood_mean <- nhood_sum %*% ncells_mat

  return(nhood_mean)

}



#' Calculates gene specificity values
calcGeneSpec <- function(nhood_mean) {
  N <- ncol(nhood_mean)
  col_sums <- Matrix::colSums(nhood_mean)
  gspec <- nhood_mean %*% diag(N/col_sums)
  return(gspec)
}



getScranHVGs <- function(sce, n_hvgs, block=NULL) {
  gene_var <- scran::modelGeneVar(sce, block=NULL)
  hvgs <- scran::getTopHVGs(gene_var,n=n_hvgs)
  return(hvgs)
}


selectNhoodFeatures <- function(milo, hvg_selection="scran", max_hvgs=2000,
                                hvg_block=NULL, exclude_genes = NULL) {

  # If hvgs already provided
  if(length(hvg_selection)>1) return(hvg_selection)

  if(!is.null(hvg_block)) {
    hvg_block <- colData(milo, hvg_block)
  }

  if(!is.null(exclude_genes)) {
    milo <- milo[!(rownames(milo) %in% exclude_genes), ]
  }

  if(hvg_selection == "scran") {
    hvgs <- getScranHVGs(milo, max_hvgs, block=hvg_block)
  }

  return(hvgs)
}


exportNhoodSim <- function(export_dir, r_vals, m_vals, nhood_sim) {
  r_vals <- as.data.frame(as.matrix(r_vals))
  m_vals = as.data.frame(as.matrix(m_vals))#

  write.table(r_vals, paste0(export_dir, "r_vals.tsv"), sep="\t", quote=FALSE)
  write.table(m_vals, paste0(export_dir, "m_vals.tsv"), sep="\t", quote=FALSE)

  write.table(nhood_sim, paste0(export_dir, "nhood_sim.tsv"), sep="\t",
              quote = FALSE)
}


#' Compare neighbourhoods across species
#'
#' Runs the neighbourhood comparison pipeline
#'
#' @param r_milo The rabbit \linkS4class{Milo} object
#' @param m_milo The mouse \linkS4class{Milo} object
#' @param orthologs A \linkS4class{DataFrame} of rabbit and mouse one-to-one orthologs.
#' @param sim_preprocessing Used to specify whether to compute the gene specificity
#'  before computing similiarities between neighbourhoods. Can be either
#'  \code{"gene_spec"} or \code{"none"}.
#' @param sim_measure The similarity measure to compare neighbourhoods. Can be
#' one of \code{"pearson"}, \code{"kendall"} or \code{"spearman"}.
#' @param hvg_join_type Specifies how to combine gene features from the two species. Can
#' be either \code{"intersection"} or \code{"union"}.
#' @param r_exclude A vector of rabbit genes to exclude from feature selection.
#' @param m_exclude A vector of mouse genes to exclude from feature selection.
#' @param export_dir A string path indicating where to export output files.
#' @param verbose A logical scalar indicating whether progress updates should
#' be printed to screen.
#' @param ... Additional arguments to pass to specific methods.
#' See \code{?selectNhoodFeatures}.
#'
#' @return
#' A named list is returned containing:
#' \describe{
#' \item{\code{r_vals}: }{A \linkS4class{DataFrame} of mean expression values (or gene specificity values) for each
#' rabbit gene across all rabbit neighbourhoods.}
#' \item{\code{m_vals}: }{A \linkS4class{DataFrame} of mean expression values (or gene specificity values) for each
#' mouse gene across all mouse neighbourhoods.}
#' \item{\code{nhood_sim}: }{An all vs all matrix of similarities between rabbit (rows) and mouse (columns) neighbourhoods.}
#' }
#'
#' @details
#' This function implements the neighbourhood comparison pipeline described in
#' Ton et al. (2022).
#'
#' The pipeline consists of two main steps - selecting features and computing neighbourhood
#' similarities.
#'
#'
#' \strong{Feature selection}
#'
#' Fistly, to remove uninformative genes from driving similarities between species,
#' a subset of genes are chosen independently for each species by calling \code{selectNhoodFeatures}.
#' By default this uses \code{scran::getScranHVGs}, however, it's also possible to
#' specify a set of predefined features by passing a list of genes to the \code{hvg_selection} argument.
#'
#' These genes are then filtered to exclude genes that are listed in \code{r_exclude} and \code{m_exclude}.
#' Only genes that are one-to-one orthologs, as specified in \code{orthologs} are retained.
#'
#' The features from each species are then combined according to \code{hvg_join_type}. This provides a
#' common gene set with which to compute similarities between neighbourhoods.
#'
#'
#' \strong{Computing neighbourhood similarities}
#'
#' Using these selected features, a mean expression profile is then computed for each neighbourhood
#' using \code{calcNhoodMean}. These expression values are extracted using the \code{logcounts} assay
#' of the rabbit and mouse \linkS4class{Milo} object which must represent normalised logcounts.
#'
#' Prior to computing the similarity between neighbourhoods an additional normalisation
#' step can be performed using the \code{sim_preprocessing = "gene_spec"} parameter option.
#' Specifically, the gene specificity (\eqn{s^i_g}{s^i_g}) is computed for each neighbourhood, given by
#' \deqn{s^i_g = \frac{g_i}{\frac{1}{N} \sum_{k=1}^N g_k}}{s^i_g = N*g_i / (g_1 + g_2 + ... + g_N)}
#' where \eqn{g_x} represents the mean expression of gene \eqn{g} in neighbourhood \eqn{x \in {1 \hdots N}}{x = 1,...,N}.
#' This can be used to account for differences in absolute values between datasets.
#'
#' Following this, the similarity between neighbourhoods is computed using the correlation
#' measure specified by \code{sim_measure}.
#'
#' The matrix of rabbit vs mouse neighbourhood similarity scores, as well as mean expression
#' values (or gene specificity values if \code{sim_preprocessing = "gene_spcec"}) are exported
#' to the directory specified by \code{export_dir}.
#'
#'
#' @author Daniel Keitley
#'
#'
#' @references
#' Ton M.L, Keitley D et al. 2022
#' Rabbit Development as a Model for Single Cell Comparative Genomics
#' \emph{Manuscript in submission}.
#'
#' @seealso \code{\link{selectNhoodFeatures}}, \code{\link{calcNhoodMean}},
#' \code{\link{calcGeneSpec}} and \code{\link{exportNhoodSim}} for more info on
#' specific steps of the pipeline.
#'
#' See an example usage \href{https://github.com/MarioniLab/RabbitGastrulation2022/blob/master/notebooks/compare_nhoods.ipynb}{here}.
#'
#' @importFrom igraph vertex_attr
#' @export
calcNhoodSim <- function(r_milo, m_milo, orthologs,
                         sim_preprocessing="gene_spec", sim_measure="pearson",
                         hvg_join_type="intersection", r_exclude = NULL, m_exclude = NULL,
			 export_dir=NULL, verbose=TRUE,
                         ...) {

  # Check data is normalised
  if(is.null(logcounts(r_milo)) | is.null(logcounts(m_milo))) {
    stop("The neighbourhood comparison pipeline requires normalised logcounts.")
  }


  # Select features
  if(verbose) message("Selecting features...")
  r_features <- selectNhoodFeatures(r_milo, exclude_genes = r_exclude, ...)
  m_features <- selectNhoodFeatures(m_milo, exclude_genes = m_exclude, ...)


  # Combine features
  if(verbose) message("Combining features...")
  if(hvg_join_type == "union") {
    sim_features <- orthologs[orthologs[,1] %in% r_features |
                                orthologs[,2] %in% m_features, ]
  } else if(hvg_join_type == "intersection") {
    sim_features <- orthologs[orthologs[,1] %in% r_features &
                                orthologs[,2] %in% m_features, ]
  } else {
    stop("Invalid hvg_join_type. Please specify either 'union' or 'intersection'.")
  }


  r_assay <- logcounts(r_milo)
  m_assay <- logcounts(m_milo)

  r_filt <- r_assay[sim_features[,1],]
  m_filt <- m_assay[sim_features[,2],]


  # Average across nhoods
  if(verbose) message("Averaging expression across neighbourhoods...")
  r_vals <- calcNhoodMean(r_filt, r_milo@nhoods)
  m_vals <- calcNhoodMean(m_filt, m_milo@nhoods)


  # Sim preprocessing
  if(!is.null(sim_preprocessing)) {
    if(sim_preprocessing == "gene_spec") {
      r_vals <- calcGeneSpec(r_vals)
      m_vals <- calcGeneSpec(m_vals)
    } else {
      stop("Unrecognised sim_processing value. Please specify one of 'gene_spec' or NULL")
    }
  }

  colnames(r_vals) <- colnames(r_milo@nhoods)
  colnames(m_vals) <- colnames(m_milo@nhoods)


  # Compute similarity
  if(verbose) message("Computing similarity across neighbourhoods...")
  nhood_sim <- matrix(0L, nrow=ncol(r_milo@nhoods), ncol(m_milo@nhoods))
  if(sim_measure %in% c("pearson", "kendall", "spearman")) {
    nhood_sim <- cor(as.matrix(r_vals), as.matrix(m_vals), method=sim_measure)
  }


  rownames(nhood_sim) <- vertex_attr(nhoodGraph(r_milo))$name
  colnames(nhood_sim) <- vertex_attr(nhoodGraph(m_milo))$name


  if(!is.null(export_dir)) {
    if(verbose) message("Exporting results...")
    exportNhoodSim(export_dir, r_vals, m_vals, nhood_sim)
  }


  return(list(r_vals = r_vals, m_vals = m_vals, nhood_sim = nhood_sim))

}


#' @importFrom igraph vertex_attr
#' @importFrom miloR nhoodGraph
#' @export
subsetMiloGraph <- function(r_milo, r_graph) {

  # Get nhood indices
  r_nhoodIDs <- as.numeric(vertex_attr(r_graph)$name)
  r_indCells <- colnames(r_milo)[r_nhoodIDs]

  # Identify cells within r_graph neighbourhoods
  rhood_hits <- Matrix::rowSums(r_milo@nhoods[,as.character(r_nhoodIDs)])

  # Milo bug workaround - index cells included in r_milo@nhoods in Milo dev branch
  rhood_cells <- unique(c(colnames(r_milo)[rhood_hits > 0], r_indCells))

  # Filter original SingleCellExperiment
  r_milo <- r_milo[,rhood_cells]

  miloR::nhoodGraph(r_milo) <- r_graph

  return(r_milo)
}


# Needed to keep track of index cells when using subseted Milo objects
# E.g. when plotting UMAP of a subsetted Milo object

#' @importFrom miloR nhoodGraph
#' @importFrom igraph vertex_attr
#' @export
addCellNamesToGraph <- function(milo) {

  nh_graph <- miloR::nhoodGraph(milo)
  nhood_ids <- as.numeric(vertex_attr(nh_graph)$name)
  nhood_inds <- colnames(milo)[nhood_ids]
  igraph::V(nh_graph)$cell_name <- nhood_inds
  miloR::nhoodGraph(milo) <- nh_graph

  return(milo)
}


#' @importFrom miloR nhoodGraph
#' @importFrom igraph vertex_attr
#' @export
addAttributeToGraph <- function(milo, id) {
  nh_graph <- miloR::nhoodGraph(milo)
  nhood_ids <- as.numeric(vertex_attr(nh_graph)$name)
  nhood_inds <- colnames(milo)[nhood_ids]
  nh_graph <- igraph::set_vertex_attr(nh_graph, id, value = colData(milo)[nhood_inds, id])
  miloR::nhoodGraph(milo) <- nh_graph
  return(milo)
}


#' @importFrom miloR nhoodGraph
#' @importFrom igraph vertex_attr
#' @export
getNhoodPositions <- function(milo, nh_graph=NULL, dimred="UMAP") {

  # Check if graph provided
  if(is.null(nh_graph)) { nh_graph <- nhoodGraph(milo) }

  # Check graph attributes
  if(is.null(vertex_attr(nh_graph)$cell_name)) {
    stop("Neighbourhood graph must have a cell_name attribute. See '?addCellNamesToGraph'. ")
  }

  # Extract nhood embedding positions
  nhoodIDs <- as.numeric(vertex_attr(nh_graph)$name)
  indCells <- vertex_attr(nh_graph)$cell_name

  nhoodPos <- data.frame(reducedDim(milo, dimred)[indCells,1:2])
  rownames(nhoodPos) <- as.character(nhoodIDs)
  colnames(nhoodPos) <- c("x","y")

  return(nhoodPos)
}




#' @export
getMaxMappings <- function(nhood_sim, nhood_axis, long_format=FALSE,
                           col.names = c("nhoods1", "nhoods2", "sim")) {

  if(!long_format) {
    nhood_sim <- reshape2::melt(as.matrix(nhood_sim))
    colnames(nhood_sim) <- col.names
  }

  dt <- as.data.table(nhood_sim)
  colnames(dt) <- col.names

  dt_key <- colnames(dt)[nhood_axis]
  max_dt <- dt[, .SD[which.max(get(col.names[3]))], by=dt_key]

  return(max_dt)
}




# Subset neighbourhoods based on colData observation
#' @importFrom miloR nhoodGraph
#' @importFrom igraph vertex_attr
subsetNhoods <- function(milo, obs, values) {
  milo_graph <- nhoodGraph(milo)
  nhood_ids <- as.numeric(vertex_attr(milo_graph)$name)
  nhood_obs <- colData(milo)[nhood_ids,obs]
  nhood_filt <- as.character(nhood_ids[nhood_obs %in% values])
  return(nhood_filt)
}


#' @importFrom igraph induced_subgraph
#' @importFrom miloR nhoodGraph
#' @export
subsetMiloGroups <- function(milo, group_by, groups) {
  milo <- addCellNamesToGraph(milo)
  milo <- addAttributeToGraph(milo, group_by)
  nhood_filt <- subsetNhoods(milo, group_by, groups)
  graph_filt <- induced_subgraph(nhoodGraph(milo), nhood_filt)
  milo_filt <- subsetMiloGraph(milo, graph_filt)

  return(milo_filt)
}


