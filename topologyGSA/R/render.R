render.significant.cliques <- function(info, alpha=0.05) {
  if (!requireNamespace("Rgraphviz", quietly=TRUE))
    stop("library Rgraphviz is missing")

  if (length(edges(info$graph)) == 0)
    stop("cannot render a graph with no edges")

  g           <- Rgraphviz::layoutGraph(info$graph)
  check       <- info$p.value < alpha
  significant <- info$cliques[check]
  pvalues     <- info$p.value[check]

  if (length(significant)) {
    nodes <- unique(unlist(significant))

    edges.rd <- matrix(data=unlist(
                         sapply(1:length(significant),
                                function(i) .clique.edges(significant[[i]], pvalues[[i]]),
                                simplify=F)),
                       ncol=2,
                       byrow=T)
    edges <- as.matrix(tapply(edges.rd[,2], edges.rd[,1], function(x) min(as.numeric(x))))

    if (min(pvalues) == 0)
      score.max <- 24
    else
      score.max <- ceiling(-log(min(pvalues)))

    if (max(pvalues) == 0)
      score.min <- 0
    else
      score.min <- floor(-log(max(pvalues)))

    score.range <- score.max-score.min
    palette     <- tim.colors(score.range+2)[2:(score.range+1)]

    colors <- apply(edges, 1, function(p) palette[ min(24, ceiling(-log(p))) - score.min ])
    names(colors) <- rownames(edges)
    edgeRenderInfo(g) <- list(col=colors)

    colors <- rep(2, length(nodes))
    names(colors) <- nodes
    nodeRenderInfo(g) <- list(fill=colors)
  }

  Rgraphviz::renderGraph(g)
  if (length(significant))
    fields::image.plot(legend.only=TRUE,
                       legend.shrink=0.3,
                       legend.args=list(text="-log(pvalue)", line=3),
                       zlim=c(score.min, score.max))
}

.clique.edges <- function(nodes, pvalue) {
  edges   <- expand.grid(nodes, nodes, stringsAsFactors=FALSE)
  diff    <- edges[,1] != edges[,2]
  edges   <- apply(edges[diff,], 1, function(r) paste(r, collapse="~"))
  pvalues <- rep(pvalue, length(edges))
  rbind(edges, pvalues)
}
