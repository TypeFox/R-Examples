cv.logit.reg.default <-
function(x,y,k,lam.vec)
{
  if (missing(k))
    stop("k cannot be missing")
  
  if (missing(lam.vec))
    stop("lam.vec cannot be missing")
    
  people <- length(y)
  parameters <- nrow(x) 
  
  error.cv <- matrix(rep(0,length(lam.vec)*k),ncol=k)
  #error.cv <- c(rep(0,length(lam.vec)))
  num.pred <- c(rep(0,length(lam.vec)))
  for (fold in 1:k){
    if (fold<k){
    	n.test <- floor(length(y)/k)
    	ind.data.test <- ((fold-1)*n.test+1):(fold*n.test)
   	} else {
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
    		out1 <- logit.reg(x,y,lam.vec[ll])
      		out <- logit.reg(data.train.x,data.train.y,lam.vec[ll])
      		if (out$nonzeros>0){
      			out2 <- logit.reg(data.train.x[out$selected,],data.train.y,0)
      			testbeta = rep(0,parameters)
				testbeta[out$selected] = out2$estimate
				est = t(data.test.x)%*%testbeta
				res = (exp(est)/(1+exp(est)) - data.test.y)^2
				sumres = sum(res)
     			error.cv[ll,fold] <- sumres
     		} else{
     			error.cv[ll,fold] = 10^10000
     		}
     		logit.out.full <- logit.reg(x,y,lambda=lam.vec[ll])
            num.pred[ll] = logit.out.full$nonzeros
    }
 }
  
  mean.error <- apply(error.cv,1,mean)
  mean.error1 <- subset(mean.error, mean.error!=Inf)
  lam.opt <- lam.vec[max(which(mean.error == min(mean.error1)))]
#   lam.opt <- lam.vec[min(which(error.cv == min(error.cv)))]

  out=list(k = k, lam.vec = lam.vec, lam.opt = lam.opt, error.cv = error.cv, mean.error=mean.error, num.pred = num.pred)
    
  class(out) <- "cv.logit.reg"
  out

}



