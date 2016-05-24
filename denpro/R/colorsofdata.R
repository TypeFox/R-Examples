colorsofdata<-function(dendat, pcf, lst, paletti=NULL, clusterlevel=NULL, nodes=NULL)
{
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs
# this version written made by Sauli Herrala

  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  i <- 1:d
  step <-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
	  
  if (is.null(paletti))
  paletti<-c("red","blue","green",
      "orange","navy","darkgreen",
      "orchid","aquamarine","turquoise",
      "pink","violet","magenta","chocolate","cyan",
      colors()[50:657],colors()[50:657])
  
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)
  # links from rec to node (invert the links in infopointer)
 
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i
  
  den2pcf<-matrix(0,n,1)
  pcf2den<-matrix(0,rnum,1)
  value<-matrix(0,n,1)
  ala <- pcf$down
  yla <- pcf$high
  alaTesti <- t(ala) * step + pcf$support[2*1:ncol(ala) -1]
  ylaTesti <- t(yla) * step + pcf$support[2*1:ncol(ala) -1]
  
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < alaTesti) | (c(dendat[i,]) >  ylaTesti)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  
  datcol<-matrix("white",n,1)
  tok <- 0
  for (i in 1:n){
      eka<-den2pcf[i]
      if (eka>0) tok<-nodefinder[eka]
      if (tok>0) datcol[i]<-col[tok]
  }
  
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}



