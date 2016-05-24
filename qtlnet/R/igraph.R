##############################################################################
##
## $Id: igraph.R,v 2012/06/06 byandell@wisc.edu Exp $
##
##     Copyright (C) 2012 Elias Chaibub Neto and Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## Routines: plot.qtlnet, graph.qtlnet, igraph.qtlnet,
##           graph.qdg, igraph.qdg, plot.qdg
##############################################################################

plot.qtlnet <- function(x, ...)
{
  gr <- igraph.qtlnet(x, ...)
  tkplot(gr, ...)
  
  invisible(gr)
}
###################################################################
graph.qtlnet <- function(x, ...) igraph.qtlnet(x, ...)
###################################################################
## This creates object of class igraph.
igraph.qtlnet <- function(x,
                         edges = get.averaged.net(x, ...),
                         loci.list = loci.qtlnet(x, ...),
                         pheno.color="green", qtl.color="red",
                         vertex.color = node.color,
                         include.qtl=TRUE,
                         ...)
{
  node.names <- levels(edges[[1]])
  if(is.null(node.names))
    node.names <- unique(c(as.character(edges[[1]]), as.character(edges[[2]])))

  if(is.null(loci.list) | !include.qtl) {
    node.color <- pheno.color
    names(edges)[3] <- "width"
  }
  else {
    loci.data.frame <- data.frame(qtl = unlist(loci.list))
    loci.data.frame$pheno <- rep(names(loci.list), sapply(loci.list, length))

    pheno.names <- node.names
    node.names <- c(pheno.names, levels(loci.data.frame[[1]]))

    edges <- cbind.data.frame(cause = c(as.character(edges[[1]]),
                                as.character(loci.data.frame[[1]])),
                              effect = c(as.character(edges[[2]]),
                                as.character(loci.data.frame[[2]])),
                              width = c(edges[[3]],
                                rep(1, nrow(loci.data.frame))))

    node.color <- rep(qtl.color, length(node.names))
    node.color[node.names %in% pheno.names] <- pheno.color
  }

  ## Not sure how these get set up and passed.
  ## Set up vertices
  vertex.color <- array(vertex.color, length(node.names))
  vertices <- data.frame(name = node.names, label = node.names,
                         color = vertex.color, fill = vertex.color)

  ## Great graph object (library igraph).
  igraph.options(print.graph.attributes = TRUE,
                 print.vertex.attributes = TRUE,
                 print.edge.attributes = TRUE)
  graph.data.frame(edges, TRUE, vertices = vertices)
}
##################################################################
## Following routines are highly dependent on how igraph objects are structured.
##################################################################
get.graph.vertices <- function(graph)
{
  attr <- list.vertex.attributes(graph)
  out <- list()
  for(i in attr)
    out[[i]] <- get.vertex.attribute(graph, i)
  data.frame(out)
}
############################################################
get.graph.edges <- function(graph)
{
  attr <- list.edge.attributes(graph)
  out <- as.data.frame(get.edgelist(graph))
  names(out) <- c("cause","effect")
  for(i in attr)
    out[[i]] <- get.edge.attribute(graph, i)
  out
}

############################################################
## Used for QDG routines.
############################################################
graph.qdg <- function(x, ...) igraph.qdg(x, ...)

############################################################
igraph.qdg <- function(x,
                       edges = myedges, loci.list = myloci.list, ...,
                       simple = FALSE)
{

  ## Prepare parameters for plotting function.
  if(inherits(x, "qdg")){ 
    best <- which(x$Solutions$BIC == min(x$Solutions$BIC))
    pheno.output <- data.frame(x$Solutions$solutions[[best]],rep(0,nrow(x$Solutions$solutions[[best]])))
  }
  else if (inherits(x, "qdg.sem")){ 
    best <- which(x$BIC.SEM[,1] == min(x$BIC.SEM[,1]))
    pheno.output <- data.frame(x$Solutions$solutions[[best]],x$path.coeffs)
  }
  names(pheno.output) <- c(names(x$Solutions$solutions[[best]]),"path")

  dir <- (pheno.output$direction == "---->")
  node1 <- as.character(pheno.output$node1)
  node2 <- as.character(pheno.output$node2)
  if(any(!dir)) {
    tmp <- node1[!dir]
    node1[!dir] <- node2[!dir]
    node2[!dir] <- tmp
  }
  loci <- x$phenotype.names
  myedges <- data.frame(cause = factor(node1, loci), effect = factor(node2, loci),
                      prob = pchisq(log(10) * pheno.output$lod, 1))
  myloci.list <- x$marker.names

  igraph.qtlnet(x, edges, loci.list, ...)
}
