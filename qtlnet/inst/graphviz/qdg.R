##############################################################################
##
## $Id: qdg.R,v 20120604 byandell@wisc.edu Exp $
##
##     Copyright (C) 2007 Elias Chaibub Neto and Brian S. Yandell
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
## Routines: graphviz.qdg, Rgraphviz.qdg
##############################################################################

## These were graph.qdg and plot.qdg in R/qdg.

graphviz.qdg <- function(x, simple = FALSE, breaks = c(1,3,10,20),
                         col = c(pos.color="green", neg.color="red", pheno.color="yellow",
                         qtl.color="magenta"),include.qtl=TRUE,
                         ...)
{
  ################################################
  create.directed.graph.object <- function(output) {
    n.edges <- length(output[,1])	
    auxDG <- data.frame(matrix(NA,n.edges,3))
    aux <- which(output[,2] == "---->")
    if(length(aux) > 0) auxDG[aux,] <- output[aux,1:3]
    aux <- which(output[,2] == "<----")
    if(length(aux) > 0){	
      auxDG[aux,1] <- output[aux,3]
      auxDG[aux,2] <- "---->"
      auxDG[aux,3] <- output[aux,1]
    }
    mynodes <- unique(c(auxDG[,1],auxDG[,3]))
    le <- length(mynodes)
    edL <- vector("list",length=le)
    names(edL) <- mynodes
    for(i in 1:le){
      auxNode <- mynodes[i]
      auxEdges <- c()
      for(j in 1:n.edges){
        if(auxDG[j,1] == auxNode){
          auxEdges <- c(auxEdges,which(mynodes==auxDG[j,3]))
        }
      }
      edL[[i]] <- list(edges = auxEdges)
    }
    new("graphNEL",nodes=mynodes,edgeL=edL,edgemode="directed")
  }
  ############################################################
  graph.and.attributes <- function(pheno.output, 
                                   qtl.output=NULL, 
                                   breaks, 
                                   pos.color="black", 
                                   neg.color="black", 
                                   pheno.color="transparent", 
                                   qtl.color="transparent", 
                                   node.shape="ellipse") {
    if(length(breaks) != 4) stop("breaks must be a vector of length 4")
    if(is.null(qtl.output)){
      gR <- create.directed.graph.object(pheno.output)
      n.edges <- length(pheno.output[,1])	
      auxDG <- data.frame(matrix(NA,n.edges,5))
      aux <- which(pheno.output[,2] == "---->")
      if(length(aux) > 0) auxDG[aux,] <- pheno.output[aux,]
      aux <- which(pheno.output[,2] == "<----")
      if(length(aux) > 0){
        auxDG[aux,1] <- pheno.output[aux,3]
        auxDG[aux,2] <- "---->"
        auxDG[aux,3] <- pheno.output[aux,1]
      }
      auxEdges <- paste(auxDG[,1],auxDG[,3],sep="~")
      auxDG1 <- data.frame(auxEdges,pheno.output[,4:5])
      auxDG1[,1] <- as.character(auxDG1[,1])
      edge.nms <- graph::edgeNames(gR)
      le <- length(edge.nms)
      auxDG <- data.frame(matrix(NA,le,3))
      names(auxDG) <- names(auxDG1)
      for(i in 1:le){
        auxDG[i,] <- auxDG1[which(auxDG1[,1] == edge.nms[i]),]
      }
      nAttrs <- list()
      eAttrs <- list()
      mynodes <- graph::nodes(gR)
      aux <- rep("black",n.edges)
      for(i in 1:n.edges){
        if(auxDG[i,3] < 0) aux[i] <- neg.color
        else if (auxDG[i,3] > 0) aux[i] <- pos.color
      } 	
      names(aux) <- auxDG[,1]
      eAttrs$color <- aux
      abs.lod <- abs(auxDG[,2])
      aux <- rep(NA,n.edges)
      aux[which(abs.lod < breaks[1])] <- 1
      aux[which(abs.lod >= breaks[1] & breaks[2] > abs.lod)] <- 2
      aux[which(abs.lod >= breaks[2] & breaks[3] > abs.lod)] <- 3
      aux[which(abs.lod >= breaks[3] & breaks[4] > abs.lod)] <- 4
      aux[which(abs.lod >= breaks[4])] <- 5
      names(aux) <- auxDG[,1]
      eAttrs$lwd <- aux
      aux <- graph::degree(gR)
      aux <- aux$inDegree+aux$outDegree
      aux <- (aux/(2*max(aux)))+0.5
      nAttrs$width <- aux
      nAttrs$height <- aux/2
      aux <- rep(pheno.color,length(mynodes))
      names(aux) <- mynodes
      nAttrs$fillcolor <- aux
    }
    else{
      gRq <- create.directed.graph.object(qtl.output)  
      output <- rbind(pheno.output,qtl.output)
      gR <- create.directed.graph.object(output)
      n.edges <- length(output[,1])	
      auxDG <- data.frame(matrix(NA,n.edges,5))
      aux <- which(output[,2] == "---->")
      if(length(aux) > 0) auxDG[aux,] <- output[aux,]
      aux <- which(output[,2] == "<----")
      if(length(aux) > 0){
        auxDG[aux,1] <- output[aux,3]
        auxDG[aux,2] <- "---->"
        auxDG[aux,3] <- output[aux,1]
      }
      auxEdges <- paste(auxDG[,1],auxDG[,3],sep="~")
      auxDG1 <- data.frame(auxEdges,output[,4:5])
      auxDG1[,1] <- as.character(auxDG1[,1])
      edge.nms <- graph::edgeNames(gR)
      le <- length(edge.nms)
      auxDG <- data.frame(matrix(NA,le,3))
      names(auxDG) <- names(auxDG1)
      for(i in 1:le){
        auxDG[i,] <- auxDG1[which(auxDG1[,1] == edge.nms[i]),]
      }
      nAttrs <- list()
      eAttrs <- list()
      mynodes <- graph::nodes(gR)
      aux <- rep("black",n.edges)
      for(i in 1:n.edges){
        if(!is.na(auxDG[i,3])){
          if(auxDG[i,3] > 0) aux[i] <- pos.color
          else if (auxDG[i,3] < 0) aux[i] <- neg.color 
        }
      }	
      names(aux) <- auxDG[,1]
      eAttrs$color <- aux
      abs.lod <- abs(auxDG[,2])
      aux <- rep(1,n.edges)
      aux[which(abs.lod < breaks[1])] <- 1
      aux[which(abs.lod >= breaks[1] & breaks[2] > abs.lod)] <- 2
      aux[which(abs.lod >= breaks[2] & breaks[3] > abs.lod)] <- 3
      aux[which(abs.lod >= breaks[3] & breaks[4] > abs.lod)] <- 4
      aux[which(abs.lod >= breaks[4])] <- 5
      names(aux) <- auxDG[,1]
      eAttrs$lwd <- aux
      aux <- graph::degree(gR)
      aux <- aux$inDegree+aux$outDegree
      auxq <- graph::degree(gRq)
      auxq <- auxq$inDegree
      auxq1 <- auxq[which(auxq > 0)]
      auxq2 <- auxq[which(auxq == 0)]
      nms <- names(aux)
      nmsq1 <- names(auxq1)
      nmsq2 <- names(auxq2)
      aux1 <- match(nmsq1,nms)
      aux2 <- match(nmsq2,nms)
      aux[aux1] <- aux[aux1]-auxq[nmsq1]
      aux[aux2] <- 0
      aux <- (aux/(2*max(aux)))+0.5
      nAttrs$width <- aux
      nAttrs$height <- aux/2
      aux <- rep(pheno.color,length(mynodes))
      aux[match(nmsq2,mynodes)] <- qtl.color
      names(aux) <- mynodes
      nAttrs$fillcolor <- aux
    }
    attrs <- list(node=list(shape=node.shape))
    mylist <- list(gR,eAttrs,nAttrs,attrs)
    names(mylist) <- c("graph","edges","nodes","all.nodes")
    mylist
  }

  ## Prepare parameters for plotting function.
  if(inherits(x, "qdgAlgo")){ 
    best <- which(x$Solutions$BIC == min(x$Solutions$BIC))
    pheno.output <- data.frame(x$Solutions$solutions[[best]],rep(0,nrow(x$Solutions$solutions[[best]])))
  }
  else if (inherits(x, "qdgSEM")){ 
    best <- which(x$BIC.SEM[,1] == min(x$BIC.SEM[,1]))
    pheno.output <- data.frame(x$Solutions$solutions[[best]],x$path.coeffs)
  }
  names(pheno.output) <- c(names(x$Solutions$solutions[[best]]),"path")
  if(simple){
    mygR <- create.directed.graph.object(pheno.output)
    mygR[[4]] <- myattrs <- list(node=list(shape="ellipse"))
  }
  else{
    if(include.qtl){
      markers <- unlist(x$marker.names)
      qtl.output <- data.frame(matrix(NA, length(markers), 5))
      qtl.output[,1] <- markers
      qtl.output[,2] <- "---->"
      qtl.output[,3] <- rep(x$phenotype.names,times=unlist(lapply(x$marker.names,length)))
      names(qtl.output) <- c("node1","direction","node2","lod","path")
      mygR <- graph.and.attributes(pheno.output = pheno.output, 
                                   qtl.output = qtl.output,
                                   breaks = breaks,
                                   pos.color = col["pos.color"],
                                   neg.color = col["neg.color"],
                                   pheno.color = col["pheno.color"],
                                   qtl.color = col["qtl.color"])
    }
    else{
      mygR <- graph.and.attributes(pheno.output = pheno.output, 
                                   qtl.output = NULL,
                                   breaks = breaks,
                                   pos.color = col["pos.color"],
                                   neg.color = col["neg.color"],
                                   pheno.color = col["pheno.color"],
                                   qtl.color = col["qtl.color"])
    }
  }
  mygR
}

Rgraphviz.qdg <- function(x, graph = graph.qdg(x, ...), ...)
{
  require(Rgraphviz)
  
  ## Plot the graph object.
  plot(mygR[[1]], edgeAttrs = mygR[[2]], nodeAttrs = mygR[[3]],
       attrs = mygR[[4]], "dot", ...)
}
