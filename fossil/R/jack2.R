`jack2` <-
function(x,taxa.row=TRUE,abund=TRUE) {
  if (taxa.row==FALSE) x<-t(x) 
  if (ncol(as.matrix(x))==1|nrow(as.matrix(x))==1) {
    x <- x
    m <- sum(x)
  }
  else if (abund==TRUE) {x <- rowSums(x) 
    m <- sum(x)}
  else {m <-length(x[1,])
    x1 <- numeric(length(x[,1]))
    for (i in 1:length(x[,1])) x1[i] <- length(x[i,][x[i,]>0])
    x<-x1}
  q1 <- length(x[x==1])
  q2 <- length(x[x==2])
  so <- length(x[x>0])
  soj <- so+q1*(2*(m-3)/m)-q2*((m-2)^2/(m*(m-1)))
  return(soj)
}

