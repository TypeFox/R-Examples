`ICE` <-
function(x,taxa.row=TRUE) {
  if (taxa.row==FALSE) x<-t(x)
  if (ncol(as.matrix(x))==1) x <- x[x>0] 
  else if (nrow(as.matrix(x))==1) x <- x[x>0]
  else {
    x1 <- numeric(length(x[,1]))
    for (i in 1:length(x[,1])) x1[i] <- length(x[i,][x[i,]>0])
    x <- x1[x1>0]
  }
  nr <- sum(x[x<=10])
  sa <- length(x[x>10])
  sr <- length(x[x<=10])
  f1 <- length(x[x==1])
  ca <- 1-(f1)/(nr)
  sumf <- 0
  for (i in 1:10) sumf <- sumf + (i*length(x[x==i]))
  g2a <- max((sr/ca)*(sumf/(nr*(nr-1)))-1,0)
  ice <- sa + sr/ca + (f1/ca)*g2a
  if (is.nan(ice)==TRUE | ice==Inf) ice<-chao2(x) 
  attr(ice,"method") <- "ICE"
  return(ice)
}

