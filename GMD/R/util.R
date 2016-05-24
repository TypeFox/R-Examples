
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Aug 25 18:08:55 EDT 2014 -0400 (Week 34)
## 
## 
## Reference: 
## 
## 
## ************************************************************************

##' Convert a BedGraph fiel to a vector of depth-like signals
##'
##' 
##' @title Convert a BedGraph fiel to a vector of depth-like signals
##' @param inFpath character, the path of a bedgraph file
##' @param chr character, the chromosome 
##' @param start numeric, the start position (0-base)
##' @param end numeric, the end position
##' @param reverse logical, whether a reverse strand
##' @return numeric
##' @author Xiaobei Zhao
bedgraph.to.depth <- function(inFpath,chr,start,end,reverse=FALSE){
  x0 <- read.table(inFpath,header=FALSE)
  x <- x0[as.character(x0[,1])==chr & x0[,2]>=start & x0[,3]<=end, ]
  v <- unlist(mapply(rep,x[,4],x[,3]-x[,2]))
  pos <- unlist(mapply(seq,x[,2]+1,x[,3]))
  ret <- rep(0,end-start)
  names(ret) <- seq(start+1,end)
  ret[as.character(pos)] <- v
  if (reverse) {
    ret <- rev(ret)
  }
  return(ret)
}


