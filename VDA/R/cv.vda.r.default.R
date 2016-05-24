cv.vda.r.default <-
function(x,y,k,lam.vec)
{
  if (missing(k))
    stop("k cannot be missing")
  
  if (missing(lam.vec))
    stop("lambda vector cannot be missing")
    
  if (length(y)!=nrow(x))
  stop("Dimention doesn't match! 
         Rows of feature matrix X must be the number of cases")
  
 y<-as.data.frame(y)
 xy<-cbind(x,y)
 xyran<-xy[sample(nrow(xy),nrow(xy)),]
 x<-xyran[,1:(ncol(xyran)-1)]
 y<-xyran[,ncol(xyran)]
  error.cv <- matrix(rep(0,length(lam.vec)*k),ncol=k)
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
    data.test.x <- x[ind.data.test,]
    data.test.y <- y[ind.data.test]
    data.train.x <- x[ind.data.train,]
    data.train.y <- y[ind.data.train]
    
    for (ll in 1:length(lam.vec)){
      vda.out.train <- vda.r(data.train.x,data.train.y,lambda=lam.vec[ll])
      class.pred.test <- predict(vda.out.train,data.test.x)
      error.cv[ll,fold] <- length(which(as.double(class.pred.test)!=data.test.y))/length(data.test.y)
    }
  }
  
  mean.error <- apply(error.cv,1,mean)
  lam.opt <- lam.vec[max(which(mean.error == min(mean.error)))]

  out=list(k = k, lam.vec = lam.vec, mean.error = mean.error, lam.opt = lam.opt, error.cv = error.cv)
    
  class(out) <- "cv.vda.r"
  out

}
