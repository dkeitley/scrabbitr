
#' Average expression values within each neighbourhood
#' @importFrom Matrix Diagonal
#' @param X Matrix of expression values
#' @param mask Indicator matrix specifying cell groupings
calcNhoodMean <- function(X, mask) {

  ncells <- colSums(mask)
  ncells_inv <- 1/ncells
  ncells_mat <- Matrix::Diagonal(x=ncells_inv)

  nhood_sum <- X %*% mask
  nhood_mean <- nhood_sum %*% ncells_mat

  return(nhood_mean)

}



#' Calculates gene specificity values
calcGeneSpec <- function(nhood_mean) {
  N <- ncol(nhood_mean)
  col_sums <- colSums(nhood_mean)
  gspec <- nhood_mean %*% diag(N/col_sums)
  return(gspec)
}



getScranHVGs <- function(sce, n_hvgs, block=NULL) {
  gene_var <- scran::modelGeneVar(sce, block=NULL)
  hvgs <- scran::getTopHVGs(gene_var,n=n_hvgs)
  return(hvgs)
}


selectNhoodFeatures <- function(milo, hvg_selection="scran", max_hvgs=2000,
                                hvg_block=NULL) {

  # If hvgs already provided
  if(length(hvg_selection)>1) return(hvg_selection)

  if(!is.null(hvg_block)) {
    hvg_block <- colData(milo, hvg_block)
  }

  if(hvg_selection == "scran") {
    hvgs <- getScranHVGs(milo, max_hvgs, block=hvg_block)
  }

  return(hvgs)
}


exportNhoodSim <- function(export_dir, r_vals, m_vals, nhood_sim) {
  r_vals <- as.data.frame(as.matrix(r_vals))
  m_vals = as.data.frame(as.matrix(m_vals))#

  write.table(r_vals, paste0(export_dir, "r_vals.tsv"), sep="\t", header=TRUE, quote=FALSE)
  write.table(m_vals, paste0(export_dir, "m_vals.tsv"), sep="\t", header=TRUE, quote=FALSE)

  write.table(nhood_sim, paste0(export_dir, "nhood_sim.tsv"), sep="\t",
              header=TRUE, quote = FALSE)
}



#' @examples
#' calNhoodSim(r_milo, m_milo, rm_orthologs, hvg_block="sample")
calcNhoodSim <- function(r_milo, m_milo, orthologs, assay="logcounts",
                         sim_preprocessing="gene_spec", sim_measure="pearson",
                         hvg_join_type="intersection", export_dir=NULL, ...) {

  # Check data is normalised
  if(is.null(logcounts(r_milo)) | is.null(logcounts(m_milo))) {
    stop("The neighbourhood comparison pipeline requires normalised logcounts.")
  }

  r_assay <- logcounts(r_milo)
  m_assay <- logcounts(m_milo)


  # Select features
  r_features <- selectNhoodFeatures(r_milo, ...)
  m_features <- selectNhoodFeatures(m_milo, ...)


  # Combine features
  if(hvg_join_type == "union") {
    sim_features <- orthologs[orthologs[,1] %in% r_features |
                                orthologs[,2] %in% m_features, ]
  } else if(hvg_join_type == "intersection") {
    sim_features <- orthologs[orthologs[,1] %in% r_features &
                                orthologs[,2] %in% m_features, ]
  } else {
    stop("Invalid hvg_join_type. Please specify either 'union' or 'intersection'.")
  }

  r_filt <- r_assay[sim_features[,1],]
  m_filt <- m_assay[sim_features[,2],]


  # Average across nhoods
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
  if(sim_measure %in% c("pearson", "kendall", "spearman")) {
    nhood_sim <- cor(as.matrix(r_vals), as.matrix(m_vals), method=sim_measure)
  }

  if(!is.null(export_dir)) {
    exportNhoodSim(export_dir, r_vals, m_vals, nhood_sim)
  }


  return(list(r_vals = r_vals, m_vals = m_vals, nhood_sim = nhood_sim))

}



subsetMiloGraph <- function(r_milo, r_graph) {

  # Get nhood indices
  r_nhoodIDs <- as.numeric(vertex_attr(r_graph)$name)
  r_indCells <- colnames(r_milo)[r_nhoodIDs]

  # Identify cells within r_graph neighbourhoods
  rhood_hits <- rowSums(r_milo@nhoods[,as.character(r_nhoodIDs)])

  # Milo bug workaround - index cells included in r_milo@nhoods in Milo dev branch
  rhood_cells <- unique(c(colnames(r_milo)[rhood_hits > 0], r_indCells))

  # Filter original SingleCellExperiment
  r_milo <- r_milo[,rhood_cells]

  return(r_milo)
}


# Needed to keep track of index cells when using subseted Milo objects
# E.g. when plotting UMAP of a subsetted Milo object
addCellNamesToGraph <- function(milo) {

  nh_graph <- nhoodGraph(milo)
  nhood_ids <- as.numeric(vertex_attr(nh_graph)$name)
  nhood_inds <- colnames(milo)[nhood_ids]
  V(nh_graph)$cell_name <- nhood_inds
  nhoodGraph(milo) <- nh_graph

  return(milo)
}



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



getMaxMappings <- function(nhood_sim, nhood_axis, long_format=FALSE) {

  if(!long_format) {
    nhood_sim <- reshape2::melt(as.matrix(nhood_sim))
    colnames(nhood_sim) <- c("nhoods1", "nhoods2", "sim")
  }

  df_sim <- nhood_sim[order(nhood_sim[,nhood_axis], -nhood_sim[,3]), ]
  max_nhoods <- df_sim[!duplicated(df_sim[,nhood_axis]),]
  rownames(max_nhoods) <- max_nhoods[,nhood_axis]
  #max_nhoods <- max_nhoods[vertex_attr(r_graph)$name,]

  return(max_nhoods)
}



subsetMiloGroups <- function(milo, group_by, groups) {

  nhood_filt <- subsetNhoods(milo, group_by, groups)
  graph_filt <- induced_subgraph(nhoodGraph(milo), nhood_filt)
  milo_filt <- subsetMiloGraph(r_milo, graph_filt)

  return(milo_filt)
}

