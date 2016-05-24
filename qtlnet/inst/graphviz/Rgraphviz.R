## Can this be rewritten to use graph.qtlnet output?
graphviz.plot.qtlnet <- function(qtlnet.object,
                                 pheno.output = summary(x)$averaged.net,
                                 marker.list = qtlnet.pheno(x, chr.pos),
                                 pheno.color="green", qtl.color="red",
                                 include.qtl=TRUE,
                                 simple = FALSE,
                                 thick=3, myfontsize=18, mynodewidth=3.0,
                                 mynodeheight1.0, plotstyle="dot",
                                 ...)
{
  if(simple){
    mygR <- create.directed.graph.object(pheno.output)
    myattrs <- list(node=list(shape="ellipse"))
    return( plot(mygR, attrs = myattrs, ...) )
  }
  else{
    if(include.qtl){
      pheno.nms <- names(marker.list)
      markers <- unlist(marker.list)
      qtl.output <- data.frame(matrix(NA, length(markers), 3))
      qtl.output[,1] <- markers
      qtl.output[,2] <- rep(pheno.nms,times=unlist(lapply(marker.list,length)))
      names(qtl.output) <- c("node1","direction","node2","prob")
      mygR <- graph.and.attributes(pheno.output = pheno.output, 
                                            qtl.output = qtl.output,
                                            pheno.color = pheno.color,
                                            qtl.color = qtl.color,
                                            thick = thick,
                                            myfontsize = myfontsize,
                                            mynodewidth = mynodewidth,
                                            mynodeheight = mynodeheight)
    }
    else{
      mygR <- graph.and.attributes(pheno.output = pheno.output, 
                                            qtl.output = NULL,
                                            pheno.color = pheno.color,
                                            qtl.color = qtl.color,
                                            thick = thick,
                                            myfontsize = myfontsize,
                                            mynodewidth = mynodewidth,
                                            mynodeheight = mynodeheight)
    }
    ## Plot the graph object.
    return( plot(mygR[[1]], edgeAttrs = mygR[[2]], nodeAttrs = mygR[[3]],
                 attrs = mygR[[4]], plotstyle) )
  }
}
################################################################
## Rgraphviz routines below.
################################################
create.directed.graph.object <- function(output)
{
  n.edges <- nrow(output)	
  auxDG <- data.frame(matrix(NA,n.edges,3))
  auxDG[,2] <- "---->"
  auxDG[,c(1,3)] <- output[,1:2]

  mynodes <- unique(c(auxDG[,1],auxDG[,3]))
  le <- length(mynodes)
  edL <- vector("list",length=le)
  names(edL) <- mynodes
  for(i in 1:le){
    auxNode <- mynodes[i]
    auxEdges <- c()
    for(j in 1:n.edges) {
      if(auxDG[j,1] == auxNode){
        auxEdges <- c(auxEdges,which(mynodes==auxDG[j,3]))
      }
    }
    edL[[i]] <- list(edges = auxEdges)
  }
  new("graphNEL",nodes=mynodes,edgeL=edL,edgemode="directed")
}
############################################################
graph.and.attributes <- function(pheno.output, qtl.output=NULL,
                                 pheno.color="transparent", qtl.color="transparent", 
                                 node.shape="ellipse", thick=3, myfontsize=1, mynodewidth, 
                                 mynodeheight)
{
  grey.scale <- grey(0:20/20)
  aux.grey <- (2/3)/20
  if(!is.null(qtl.output)){
    qtl.output[,4] <- 1
    pheno.output <- rbind(qtl.output,pheno.output) 
  }
  gR <- create.directed.graph.object(pheno.output)
  n.edges <- length(pheno.output[,1])	
  auxDG <- data.frame(matrix(NA,n.edges,4))
  aux <- which(pheno.output[,2] == "---->")
  if(length(aux) > 0) auxDG[aux,] <- pheno.output[aux,]
  aux <- which(pheno.output[,2] == "<----")
  if(length(aux) > 0){
    auxDG[aux,1] <- pheno.output[aux,3]
    auxDG[aux,2] <- "---->"
    auxDG[aux,3] <- pheno.output[aux,1]
    auxDG[aux,4] <- pheno.output[aux,4]
  }
  auxEdges <- paste(auxDG[,1],auxDG[,3],sep="~")
  auxDG1 <- data.frame(auxEdges,pheno.output[,4])
  auxDG1[,1] <- as.character(auxDG1[,1])
  edge.nms <- edgeNames(gR)
  le <- length(edge.nms)
  
  nAttrs <- list()
  eAttrs <- list()
  mynodes <- nodes(gR)
  aux <- rep(NA,n.edges)
  for(i in 1:n.edges){
    if((pheno.output[i,4] >= 1/3) & (pheno.output[i,4] < 1/3+1*aux.grey)) 
      aux[i] <- grey.scale[20]
    if((pheno.output[i,4] >= 1/3+1*aux.grey) & (pheno.output[i,4] < 1/3+2*aux.grey)) 
      aux[i] <- grey.scale[19]
    if((pheno.output[i,4] >= 1/3+2*aux.grey) & (pheno.output[i,4] < 1/3+3*aux.grey)) 
      aux[i] <- grey.scale[18]
    if((pheno.output[i,4] >= 1/3+3*aux.grey) & (pheno.output[i,4] < 1/3+4*aux.grey)) 
      aux[i] <- grey.scale[17]
    if((pheno.output[i,4] >= 1/3+4*aux.grey) & (pheno.output[i,4] < 1/3+5*aux.grey)) 
      aux[i] <- grey.scale[16]
    if((pheno.output[i,4] >= 1/3+5*aux.grey) & (pheno.output[i,4] < 1/3+6*aux.grey)) 
      aux[i] <- grey.scale[15]
    if((pheno.output[i,4] >= 1/3+6*aux.grey) & (pheno.output[i,4] < 1/3+7*aux.grey)) 
      aux[i] <- grey.scale[14]
    if((pheno.output[i,4] >= 1/3+7*aux.grey) & (pheno.output[i,4] < 1/3+8*aux.grey)) 
      aux[i] <- grey.scale[13]
    if((pheno.output[i,4] >= 1/3+8*aux.grey) & (pheno.output[i,4] < 1/3+9*aux.grey)) 
      aux[i] <- grey.scale[12]
    if((pheno.output[i,4] >= 1/3+9*aux.grey) & (pheno.output[i,4] <= 1/3+10*aux.grey)) 
      aux[i] <- grey.scale[11]
    if((pheno.output[i,4] >= 1/3+10*aux.grey) & (pheno.output[i,4] < 1/3+11*aux.grey)) 
      aux[i] <- grey.scale[10]
    if((pheno.output[i,4] >= 1/3+11*aux.grey) & (pheno.output[i,4] < 1/3+12*aux.grey)) 
      aux[i] <- grey.scale[9]
    if((pheno.output[i,4] >= 1/3+12*aux.grey) & (pheno.output[i,4] < 1/3+13*aux.grey)) 
      aux[i] <- grey.scale[8]
    if((pheno.output[i,4] >= 1/3+13*aux.grey) & (pheno.output[i,4] < 1/3+14*aux.grey)) 
      aux[i] <- grey.scale[7]
    if((pheno.output[i,4] >= 1/3+14*aux.grey) & (pheno.output[i,4] < 1/3+15*aux.grey)) 
      aux[i] <- grey.scale[6]
    if((pheno.output[i,4] >= 1/3+15*aux.grey) & (pheno.output[i,4] < 1/3+16*aux.grey)) 
      aux[i] <- grey.scale[5]
    if((pheno.output[i,4] >= 1/3+16*aux.grey) & (pheno.output[i,4] < 1/3+17*aux.grey)) 
      aux[i] <- grey.scale[4]
    if((pheno.output[i,4] >= 1/3+17*aux.grey) & (pheno.output[i,4] < 1/3+18*aux.grey)) 
      aux[i] <- grey.scale[3]
    if((pheno.output[i,4] >= 1/3+18*aux.grey) & (pheno.output[i,4] < 1/3+19*aux.grey)) 
      aux[i] <- grey.scale[2]
    if((pheno.output[i,4] >= 1/3+19*aux.grey) & (pheno.output[i,4] <= 1.00)) 
      aux[i] <- grey.scale[1]
  }
  names(aux) <- as.character(auxDG1[,1])
  eAttrs$color <- aux
  aux <- rep(thick,n.edges)
  names(aux) <- as.character(auxDG1[,1])
  eAttrs$lwd <- aux
  le <- length(unique(qtl.output[,1]))
  aux <- c(rep(qtl.color,le),rep(pheno.color,length(mynodes)-le))
  names(aux) <- as.character(mynodes)
  nAttrs$fillcolor <- aux
  aux <- rep(myfontsize,length(mynodes))
  names(aux) <- as.character(mynodes)
  nAttrs$fontsize <- aux
  aux <- rep(mynodewidth,length(mynodes))
  names(aux) <- as.character(mynodes)
  nAttrs$width <- aux
  aux <- rep(mynodeheight,length(mynodes))
  names(aux) <- as.character(mynodes)
  nAttrs$height <- aux
      attrs <- list(node=list(shape=node.shape))
  mylist <- list(gR,eAttrs,nAttrs,attrs)
  names(mylist) <- c("graph","edges","nodes","all.nodes")
  mylist
}
