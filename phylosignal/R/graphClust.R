

#' Phylogenetically constrained  clustering
#'
#' This function extracts clusters of species based on traits values and phylogenetic proximities.
#' 
#' @param p4d a \code{phylo4d} object.
#' @param trait the traits in the \code{phylo4d} object to use for clustering.
#' Can be a character vector giving the name of the traits or numbers giving the column index
#' in the table of the data slot of the \code{phylo4d} object.
#' @param lim.phylo the maximum phylogenetic distance for edges selection.
#' @param lim.trait the maximum trait-based distance for edges selection.
#' @param dist.phylo a matrix of phylogenetic distances or a character string specifying a method to compute it.
#' See Details.
#' @param dist.trait a character string specifying the method used to compute traits distances.
#' See Details.
#' @param select.method a character string specifying the method used to select edges.
#' This must be one of "\code{line}", "\code{rectangle}" or "\code{ellipse}".
#' @param scale.lim logical (default \code{TRUE}) indicating if \code{lim.phylo} and \code{lim.trait} are scaled
#' (divided by their max value).
#' 
#' @details If "\code{dist.phylo}" is a character string,
#' the phylogenetic distance matrix is computed internally
#' using the function \code{\link[adephylo]{distTips}} from the package \pkg{adephylo}.
#' All the methods supported by \code{\link[adephylo]{distTips}} are available:
#' "\code{patristic}","\code{nNodes}","\code{Abouheif}" and "\code{sumDD}".
#' See \code{\link[adephylo]{distTips}} for details about the methods.
#' 
#' If "\code{dist.trait}" is a character string,
#' the traits distance matrix is computed with the \code{\link[stats]{dist}} function.
#' All the methods supported by \code{\link[stats]{dist}} are available:
#' "\code{euclidean}","\code{maximum}", "\code{manhattan}",
#' "\code{canberra}", "\code{binary}" and "\code{minkowski}".
#' See \code{\link[stats]{dist}} for details about the methods.
#' 
#' @return An object of class \code{graphclust}.
#'
#' @examples
#' data(navic)
#' gC <- graphClust(navic, lim.phylo = 1, lim.trait = 2, scale.lim = FALSE)
#' gC
#' plot.graphclust(gC, which = "selection", ask = FALSE)
#' plot.graphclust(gC, which = "graph", ask = FALSE)
#' plot.graphclust(gC, which = "tree", ask = FALSE)
#' 
#' @export
graphClust <- function(p4d, trait = names(tdata(p4d)),
                       lim.phylo = 0.2, lim.trait = 0.2, select.method = "ellipse",
                       dist.phylo = "patristic", dist.trait = "euclidean",
                       scale.lim = TRUE){
  
  select.method <- match.arg(select.method, c("line", "rectangle", "ellipse"))
  
  if(is.numeric(trait)){
    trait <- names(tdata(p4d))[trait]
  }
  
  p4 <- extractTree(p4d)
  phy <- as(p4, "phylo")
  new.order <- phy$edge[, 2][!phy$edge[, 2] %in% phy$edge[, 1]]
  tips <- phy$tip.label[new.order]
  n.tips <- length(tips)
  X <- tdata(p4d, type = "tip")
  X <- X[tips, trait]
  X <- scale(X)
  X <- as.data.frame(X)
  colnames(X) <- trait
  rownames(X) <- tips
  n.traits <- ncol(X)
  
  if(is.vector(dist.phylo) & is.character(dist.phylo)){
    dist.phylo <- match.arg(dist.phylo, c("patristic", "nNodes", "Abouheif", "sumDD"))
    dmat.phylo <- distTips(phy, method = dist.phylo)
    dmat.phylo <- as.matrix(dmat.phylo)
    dmat.phylo <- dmat.phylo[tips, tips]
  } else {
    if(is.matrix(dist.phylo)){
      dmat.phylo <- dist.phylo[tips, tips]
    } else {
      stop("dist.phylo is not valid")
    }
  }
  
  if(is.vector(dist.trait) & is.character(dist.trait)){
    dmat.trait <- dist(X, method = dist.trait)
    dmat.trait <- as.matrix(dmat.trait)
    dmat.trait <- dmat.trait[tips, tips]
  } else {
    if(is.matrix(dist.trait)){
      dmat.trait <- dist.trait[tips, tips]
    } else {
      stop("dist.trait is not valid")
    }
  }
  
  if(scale.lim){
    lim.p <- max(dmat.phylo) * lim.phylo
    lim.t <- max(dmat.trait) * lim.trait
  } else {
    lim.p <- lim.phylo
    lim.t <- lim.trait
  }
  
  if(select.method == "line"){
    adj.mat <- (((-lim.t/lim.p) * dmat.phylo) + lim.t) > dmat.trait
  }
  if(select.method == "rectangle"){
    adj.mat <- (dmat.phylo < lim.p) + (dmat.trait < lim.t) == 2
  }
  if(select.method == "ellipse"){
    adj.mat <- (dmat.phylo/lim.p)^2 + (dmat.trait/lim.t)^2 < 1
  }
  
  gr <- igraph::graph.adjacency(adj.mat, mode = "lower", diag = FALSE)
  
  gr.decomp <- igraph::decompose.graph(gr)
  gr.decomp.names <- lapply(gr.decomp, function(x) igraph::V(x)$name)
  gr.decomp.den <- sapply(gr.decomp, igraph::graph.density)
  names(gr.decomp.den) <- 1:length(gr.decomp)
    
  clust.gr <- igraph::clusters(gr)
  clust <- clust.gr$membership
  names(clust) <- igraph::V(gr)$name
  
  inclusion <- (igraph::degree(gr) + 1) / clust.gr$csize[clust.gr$membership]
  
  meta <- list(p4d = p4d, trait = trait, adj.mat = adj.mat,
               dmat.phylo = dmat.phylo, dmat.trait = dmat.trait,
               lim.p = lim.p, lim.t = lim.t,
               select.method = select.method,
               dist.phylo = dist.phylo, dist.trait = dist.trait,
               graph = gr, clust.gr = clust.gr)
  res <- list(clusters = clust,
              clusters.density = gr.decomp.den,
              taxa.inclusion = inclusion,
              meta = meta)
  
  class(res) <- "graphclust"
  return(res)
}



#' Plot phylogenetically constrained  clustering
#'
#' This function produces three plots (selectable by \code{which}):
#' a plot of edges selection based on phylogenetic against trait distances of taxa pairs,
#' a plot of the graph produced with the selected edges
#' and a plot of the clustered phylogenetic tree.
#' 
#' @param x a \code{graphclust} object as produced by \code{\link{graphClust}}.
#' @param which a character vector to select plots.
#' Must be one or more of "\code{selection}", "\code{graph}", "\code{tree}".
#' @param ask logical if \code{TRUE} (default), the user is asked before each plot.
#' @param colored logical indicating if plots must include colors.
#' @param ... further arguments to be passed to or from other methods. They are currently ignored in this function.
#'
#' @examples
#' data(navic)
#' gC <- graphClust(navic, lim.phylo = 1, lim.trait = 2, scale.lim = FALSE)
#' plot.graphclust(gC, which = "selection", ask = FALSE)
#' plot.graphclust(gC, which = "graph", ask = FALSE)
#' plot.graphclust(gC, which = "tree", ask = FALSE)
#'
#' @method plot graphclust
#' @export
plot.graphclust <- function(x, which = c("selection", "graph", "tree"), ask = TRUE, colored = TRUE, ...){
  if(class(x) != "graphclust"){
    stop("x must be an object of class 'graphclust'")
  }
  which <- match.arg(which, c("selection", "graph", "tree"), several.ok = TRUE)
  ask0 <- par("ask")
  
  if(x$meta$select.method == "line"){
    gx <- c(0, x$meta$lim.p)
    gy <- c(x$meta$lim.t, 0)
  }
  if(x$meta$select.method == "rectangle"){
    gx <- c(0, x$meta$lim.p, x$meta$lim.p)
    gy <- c(x$meta$lim.t, x$meta$lim.t, 0)
  }
  if(x$meta$select.method == "ellipse"){
    gx <- seq(0, x$meta$lim.p, length.out = 100)
    gy <- sqrt((1-((gx^2)/(x$meta$lim.p^2)))*(x$meta$lim.t^2))
  }
  
  if(ask) par(ask = TRUE)
  if("selection" %in% which){
    if(colored){
      dot.col <- x$meta$adj.mat + 1
      line.col <- 2
    } else {
      dot.col <- 1
      line.col <- 1
    }
  plot(as.vector(x$meta$dmat.phylo), as.vector(x$meta$dmat.trait),
       col = dot.col,
       pch = ifelse(x$meta$adj.mat == 0, 1, 4),
       main = "Selected Edges",
       xlab = paste("Phylogenetic distance (", x$meta$dist.phylo, ")", sep = ""),
       ylab = paste("Trait distance (", x$meta$dist.trait, ")", sep = ""))
  mtext(paste0(sum(x$meta$adj.mat[lower.tri(x$meta$adj.mat, diag = FALSE)]), "/",
               length(x$taxa.inclusion) * length(x$taxa.inclusion) - length(x$taxa.inclusion)),
        line = 0.5, side= 3)
  lines(gx, gy, col = line.col, lwd = 2, lty = "dashed")
  }
  if(colored){
    graph.gcol <- tree.gcol <- evenColors(x$meta$clust.gr$no)[x$clusters]
    names(graph.gcol) <- names(tree.gcol) <- names(x$clusters)
  } else {
    graph.gcol <- "grey"
    tree.gcol <- "grey35"
  }
  
  if("graph" %in% which){
    igraph::plot.igraph(x$meta$graph, vertex.color = graph.gcol, vertex.label.color = 1, vertex.label.family = "")
  }
  if("tree" %in% which){
    barplot(x$meta$p4d, x$meta$trait, bar.col = tree.gcol, center = T, scale = T)
  }
  par(ask = ask0)
}

#' #' Print results of phylogenetically constrained  clustering
#' 
#' @param x a \code{graphclust} object as produced by \code{\link{graphClust}}.
#' @param ... further arguments to be passed to or from other methods.
#' 
#' @method print graphclust
#' @export
print.graphclust <- function(x, ...){
  cat("Phylogenetically constrained clustering of", length(x$taxa.inclusion), "taxa\n")
  cat("\t Included trait(s):", paste(x$meta$trait, collapse=", "), "\n")
  cat("\t Phylogenetic distance:", x$meta$dist.phylo, "\n")
  cat("\t Trait distance:", x$meta$dist.trait, "\n")
  cat("\t Selected edges:", sum(x$meta$adj.mat[lower.tri(x$meta$adj.mat, diag = FALSE)]), "\n")
  cat("\t Number of clusters:", x$meta$clust.gr$no)
}

