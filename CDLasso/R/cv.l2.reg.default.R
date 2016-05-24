cv.l2.reg.default <-
function(x,y,k,lam.vec)
{
  if (missing(k))
    stop("k cannot be missing")
  
  if (missing(lam.vec))
    stop("lam.vec cannot be missing")

  error.cv <- matrix(rep(0,length(lam.vec)*k),ncol=k)
  num.pred <- rep(0,length(lam.vec))
  
  for (fold in 1:k){
    if (fold<k){
      n.test <- floor(length(y)/k)
      ind.data.test <- ((fold-1)*n.test+1):(fold*n.test)
    }else{
      n.test <- length(y)-floor(length(y)/k)*(k-1)
      ind.data.test <- (length(y)-n.test+1):(length(y))
    }
    n.train <- length(y)-n.test
    ind.data.train <- which((1:length(y)) %in% ind.data.test==FALSE)
    data.test.x <- x[,ind.data.test]
    data.test.y <- y[ind.data.test]
    data.train.x <- x[,ind.data.train]
    data.train.y <- y[ind.data.train]
    
    for (ll in 1:length(lam.vec)){
    	eps = 0.001
      L2.out.train <- l2.reg(data.train.x,data.train.y,lambda=lam.vec[ll])
      #num.pred[ll] = L2.out.train$nonzeros
      L2.out.test <- sum(abs(t(data.test.x)%*%L2.out.train$estimate-data.test.y))
      error.cv[ll,fold] <- L2.out.test
      L2.out.full <- l2.reg(x,y,lambda=lam.vec[ll])
      num.pred[ll] = L2.out.full$nonzeros

    }
  }
  
  mean.error <- apply(error.cv,1,mean)
  lam.opt <- lam.vec[max(which(mean.error == min(mean.error)))]

  out=list(k = k, lam.vec = lam.vec, mean.error = mean.error, lam.opt = lam.opt, error.cv = error.cv, num.pred=num.pred)
    
  class(out) <- "cv.l2.reg"
  out

}
