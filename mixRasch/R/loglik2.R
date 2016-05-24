`loglik2` <-
function(Pxji, x.i, ignore.extreme = FALSE, i.stat, treat.extreme){
  n.i  <- ncol(x.i) 
  n.unique <- nrow(x.i)
  x.i <- x.i + 1

  llik <- rep(0,n.unique)

  for(i in 1:n.i){
    if(ignore.extreme){
     if(i.stat$Si[i] > treat.extreme & (i.stat$n.ni[i] - i.stat$Si[i]) > treat.extreme){
      Pxj   <- matrix(Pxji[,,i],nrow=dim(Pxji)[1])
      all.p <- rbind(1-colSums(Pxj,na.rm=TRUE),Pxj)
      all.p <- log(all.p[cbind(x.i[,i],1:n.unique)])
      llik <- llik + ifelse(is.na(all.p),0,all.p)
     }
    } else{
    Pxj   <- matrix(Pxji[,,i],nrow=dim(Pxji)[1])
    all.p <- rbind(1-colSums(Pxj,na.rm=TRUE),Pxj)
    all.p <- log(all.p[cbind(x.i[,i],1:n.unique)])
    llik <- llik + ifelse(is.na(all.p),0,all.p)
    }
  }

  llik
}

