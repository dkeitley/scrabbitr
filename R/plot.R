

#'
#'@export
plotPseudocells <- function(sce, pseudo_embed, dimred="UMAP",
                            colour_by="celltype",colours=celltype_colours,
                            point_size=1) {

  df <- reducedDim(sce,dimred)
  df <- as.data.frame(df)
  df <- cbind(df, sce[[colour_by]])
  colnames(df) <- c("x","y",colour_by)

  pseudo_embed <- as.data.frame(pseudo_embed)
  colnames(pseudo_embed) <- c("x","y")


  p <- ggplot2::ggplot(df, aes_string(x="x",y="y",colour=colour_by),
               size=point_size) +
    scale_color_manual(values = colours, name = "") +
    ggrastr::geom_point_rast(alpha=0.8,stroke=0,shape=16) +
    geom_point(data=pseudo_embed,aes(x=x,y=y),color="black",
               fill='red',stroke =0.5,size=point_size*1.5,shape=21) +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    xlab(paste0(dimred,"_1")) + ylab(paste0(dimred,"_2")) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          aspect.ratio = 1)


  return(p)

}





#' @importFrom stats aggregate
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggalluvial geom_alluvium geom_stratum geom_flow
#' @importFrom ggplot2 ggplot geom_text scale_x_discrete theme
#' @importFrom ggplot2 scale_fill_manual theme_minimal aes element_blank
#' @importFrom ggplot2 element_text
#' @export
plotAnnotationAlluvium <- function(sce, old, new, ncell_thresh=10,
                                   repel_labels=TRUE, stratum_labels=FALSE,
                                   nudge=0.2, palette = getCelltypeColours()) {

  df_frac <- aggregate(x = colnames(sce), by = list(sce[[old]],sce[[new]]),
                       FUN = length)
  colnames(df_frac) <- c("original", "predicted", "num_cells")
  df_frac[["original"]] <- as.factor(paste0("orig_",df_frac[["original"]]))
  df_frac[["predicted"]] <- as.factor(paste0("pred_", df_frac[["predicted"]]))


  # Filter entries according to the number of cells
  df_frac <- df_frac[df_frac$num_cells > ncell_thresh,]

  names(palette) <- paste0("pred_",names(palette))


  p <- ggplot(df_frac,aes(axis1 = original, axis2 = predicted, y=num_cells, fill=predicted, label=predicted)) +
    geom_alluvium() +  geom_stratum() +
    scale_x_discrete(limits = c("Original annotation", "Predicted annotation"),
                     expand = c(.4, 0.4)) +
    geom_flow()

  if(repel_labels) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = ifelse(substr(as.character(after_stat(stratum)), start = 1, stop = 5) == "orig_",
                         substring(as.character(after_stat(stratum)), 6) , NA)),
      stat = "stratum", size = 4, direction = "y", nudge_x = -nudge
    ) +

      ggrepel::geom_text_repel(
        aes(label = ifelse(substr(as.character(after_stat(stratum)), start = 1, stop = 5)=="pred_",
                           substring(as.character(after_stat(stratum)), 6), NA)),
        stat = "stratum", size = 4, direction = "y", nudge_x = nudge
      )
  }

  if(stratum_labels) {
    p <- p + geom_text(stat = "stratum",
                       aes(label = substring(as.character(after_stat(stratum)), 6)),
                       min.y=30)
  }

  p <- p + theme_minimal() +
    ylab("Number of cells") +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.background = element_blank(),
      text = element_text(size=20)
    ) +
    scale_fill_manual(values=palette) +
    theme(legend.position = "none")

  return(p)

}




#' @importFrom miloR nhoodGraph
#' @importFrom igraph simplify vertex_attr V
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggraph ggraph geom_node_point geom_edge_link0
#' @importFrom ggplot2 theme aes theme_classic element_blank
#' @importFrom viridis scale_color_viridis
#' @importFrom ggrastr rasterise
#' @export
plotNhoodMaxSim <- function(milo, df_maxNhood, colour_by="sim", legend_title="Max correlation" ) {

  nh_graph <- nhoodGraph(milo)
  V(nh_graph)$max_correlation <- df_maxNhood[vertex_attr(nh_graph)$name, colour_by]
  layout <- reducedDim(milo, "UMAP")[as.numeric(vertex_attr(nhoodGraph(milo))$name),]
  colnames(layout) <- c("x","y")

  p <- ggraph(simplify(nh_graph), layout = layout) +
    rasterise(geom_edge_link0(edge_colour = "grey66", edge_alpha=0.1),dpi=300) +
    rasterise(geom_node_point(aes(color = max_correlation), size = 3,alpha=0.8, shape=20),dpi=300) +
    theme_classic(base_size=14) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          aspect.ratio=1) +
    scale_color_viridis(option="cividis",name="Max correlation")

  return(p)

}



#' sim_values in order of milo nhoods
#' TODO: Make into separate functions so can pass plot specific params using ...
#'
#' @importFrom miloR nhoodGraph
#' @importFrom SingleCellExperiment colData
#' @importFrom igraph vertex_attr
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin aes coord_flip facet_grid
#' @importFrom ggplot2 theme guides scale_fill_manual xlab ylab
#' @importFrom ggridges geom_density_ridges
#' @export
plotNhoodSimGroups <- function(milo, sim_values, group_by="celltype",facet_by=NULL, type="ridge",
                               orientation="vertical",decreasing=FALSE,rel_min_height=0,
                               group_colours=celltype_colours,
                               xlabel="Cell type", ylabel="Correlation") {
  graph <- nhoodGraph(milo)

  df <- data.frame(nhood=vertex_attr(graph)$name,
                   max_sim=sim_values,
                   group=colData(milo)[as.numeric(vertex_attr(graph)$name), group_by])

  if(!is.null(facet_by)) {
    df$facet_by <- colData(milo)[as.numeric(vertex_attr(graph)$name), facet_by]
  }

  df_mean <- aggregate(df[,"max_sim"], list(df$group), mean,)
  df$group <- factor(df$group, levels=df_mean[order(df_mean$x,decreasing=decreasing),"Group.1"])

  p <- ggplot(df, aes(x=max_sim,y=group))

  if(type=="boxplot") {
    p <- p + geom_boxplot(aes(fill=group),outlier.size=0.1)

  } else if(type=="violin") {
    p <- p + geom_violin(aes(fill=group))

  } else {
    p <- p + geom_density_ridges(aes(fill = group),size=0.15,rel_min_height = rel_min_height)
  }


  if(!is.null(facet_by)) p <- p + facet_grid(cols=vars(facet_by))

  if(orientation=="horizontal") {
    p <- p + coord_flip()
    if(!is.null(facet_by)) p <- p + facet_grid(rows=vars(facet_by))
  }

  p <- p + theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill="none") +
    scale_fill_manual(values = group_colours[names(group_colours) %in% unique(df$group)], name = "") +
    xlab(xlabel) + ylab(ylabel)

  return(p)

}



.calcEuclidean <- function(x1,y1,x2,y2) (x1-y1)**2 + (x2-y2)**2


#' Plot a UMAP with annotated labels
#' Used to show predicted cell type annotations, where the number of factors is
#' large.
#'
#' @export
plotAnnotatedUMAP <- function(sce, colour_by="predicted_celltype", ncell_filt=5,
                              size=1, label_force=15, text_size = 2, line_size = 0.5,
                              palette=NULL) {

  if(is.null(palette)) palette <- getCelltypeColours()

  col_data <- as.data.frame(colData(sce))
  col_data$group <- col_data[[colour_by]]

  umap_df <- as.data.frame(reducedDim(sce,"UMAP"))
  colnames(umap_df) <- c("UMAP_1","UMAP_2")
  col_data <- cbind(col_data, umap_df)

  # Get mean position for each group
  mean_data <- col_data %>% group_by(group) %>% summarize_at(.vars = vars(UMAP_1,UMAP_2),.funs = c("mean"))
  mean_data <- as.data.frame(mean_data[complete.cases(mean_data),])
  rownames(mean_data) <- mean_data$group

  # Get position of closest cell to group mean
  label_pos <- col_data %>% group_by(group) %>%  filter(
    .calcEuclidean(UMAP_1, mean_data[group,"UMAP_1"], UMAP_2, mean_data[group,"UMAP_2"]) ==
      min(.calcEuclidean(UMAP_1, mean_data[group,"UMAP_1"], UMAP_2, mean_data[group,"UMAP_2"])))

  # Filter annotations with less than ncell_filt cells
  freqs <- table(sce[[colour_by]])
  label_pos <- label_pos[label_pos$group %in% names(freqs[freqs > ncell_filt]),]

  # Repels labels from rest of points
  col_data$group <- ""
  label_pos <- rbind(label_pos, col_data)

  # Wrap long labels
  label_pos$group_wrapped <- stringr::str_wrap(label_pos$group , width = 10)

  p <- ggplot(col_data, aes_string(x="UMAP_1",y="UMAP_2",colour=colour_by)) +
    ggrastr::geom_point_rast(size=size) +
    geom_text_repel(data=label_pos, aes(x=UMAP_1, y=UMAP_2,label=group_wrapped,segment.colour=group),color="black",
                    min.segment.length = 0,box.padding = 0.5,max.overlaps=Inf,size=text_size,force=label_force,
		    segment.size = line_size) +
    coord_cartesian(clip = "off") +
    scale_colour_manual(aesthetics=c("color","segment.colour"),values=palette[sce[[colour_by]]],drop=TRUE,
                        breaks=names(freqs[freqs > ncell_filt])) +

    theme_void()+
    theme(legend.position="none")

  return(p)


}

#' @importFrom viridis scale_colour_viridis
#' @export
plotGeneUMAP <- function(sce, gene, point_shape=16, point_size=1, opacity=0.8,
                         rasterise = TRUE) {
  p <- ggcells(sce, aes_string(x="UMAP.1",y="UMAP.2",colour=gene))

  if(rasterise) p <- p + ggrastr::geom_point_rast(size=point_size, alpha=opacity,
                                                  shape=point_shape)
  else p <- p + geom_point(size=point_size,alpha=opacity,shape=point_shape)

  p <- p + scale_colour_viridis(direction=-1, name = gene) +
    xlab("UMAP1") + ylab("UMAP2") +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.text = element_blank(),
          axis.ticks=element_blank(),
          aspect.ratio=1)

  return(p)
}


