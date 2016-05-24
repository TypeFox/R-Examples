colorsofdata.adagrid <- function(dendat, pcf, lst, paletti = NULL, clusterlevel=NULL, nodes=NULL){
  # links from dendat to rec to node to color
  # "lst$infopointer" gives links from nodes to recs
# this version written made by Sauli Herrala

  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  	  
  if (is.null(paletti)){
    paletti<-c("red","blue","green","orange","navy","darkgreen",
      "orchid","aquamarine","turquoise", "pink","violet","magenta",
      "chocolate","cyan", colors()[50:657],colors()[50:657])
  }
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)

  # links from rec to node (invert the links in infopointer)
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

  # find links from dendat to rec
  den2pcf<-matrix(0, n, 1)
  pcf2den<-matrix(0, rnum, 1)
  value<-matrix(0, n, 1)
  ala <- pcf$down
  yla <- pcf$high
  # a bit complex
  ala <- sapply(1:ncol(ala), function(x, pcf, ala) pcf$grid[ala[, x], x], pcf = pcf, ala = ala)
  yla <- sapply(1:ncol(yla), function(x, pcf, yla) pcf$grid[yla[, x], x], pcf = pcf, yla = yla)
  ala <- t(ala)
  yla <- t(yla)
  
  prc <- proc.time()
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < ala) | (c(dendat[i,]) >  yla)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  datcol <- matrix("white",n,1)
  for (i in 1:n){
    eka<-den2pcf[i]
    if (eka>0) tok<-nodefinder[eka]
    if ((eka>0)&&(tok>0)) datcol[i]<-col[tok]
  } 
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}

