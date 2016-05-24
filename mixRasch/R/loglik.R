`loglik` <-
function(Pxji, x.i){
  n.i  <- dim(Pxji)[3] 
  n.unique <- nrow(x.i)
  x.i <- x.i + 1

  llik <- rep(0,n.unique)

  for(i in 1:n.i){
      Pxj   <- matrix(Pxji[,,i],nrow=dim(Pxji)[1])
      all.p <- rbind(1-colSums(Pxj,na.rm=TRUE),Pxj)
      all.p <- log(all.p[cbind(x.i[,i],1:n.unique)])
      llik <- llik + ifelse(is.na(all.p),0,all.p) 
  }

  llik
}

