`ACE` <-
function(x,taxa.row=TRUE) {
  if (taxa.row==FALSE) x<-t(x)
  if (ncol(as.matrix(x))==1 | nrow(as.matrix(x))==1) x <- x[x>0] 
  else {x <- rowSums(x)
  x <- x[x>0]}
  nr <- sum(x[x<=10])
  sa <- length(x[x>10])
  sr <- length(x[x<=10])
  f1 <- length(x[x==1])
  ca <- 1-(f1)/(nr)
  sumf <- 0
  for (i in 1:10) sumf <- sumf + (i*length(x[x==i]))
  g2a <- max((sr/ca)*(sumf/(nr*(nr-1)))-1,0)
  ace <- sa + sr/ca + (f1/ca)*g2a
  if (is.nan(ace)==TRUE | ace==Inf) ace<-chao1(x) 
  attr(ace,"method") <- "ACE"
  if (sum(x[x>1])==0) warning("This data appears to be presence/absence based, but this estimator is for abundance data only")
  return(ace)
}

