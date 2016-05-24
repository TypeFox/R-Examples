predict.vda.r <-
function(object, newdata = NULL, ...)
{
  
  if (!inherits(object, "vda.r")) 
    stop("predict.VDA can only be used to predict from VDA objects")
  
  if (missing(newdata))
    return(object$predicted)
  else {
  	if (nrow(newdata)==1)  {
  		newdata <- as.matrix(newdata,nrow=1)
  		feature <- as.matrix(object$feature[,2:ncol(object$feature)])
    	newdata2 <- rbind(newdata,feature) 
    	newdata2 <- apply(newdata2,2,scale)
 
		pred <- cbind(rep(1,nrow(newdata2)),newdata2)%*%t(object$coefficient)		
 		}
 	if (nrow(newdata)>1) {
 		newdata <- as.matrix(newdata)
  		feature <- as.matrix(object$feature[,2:ncol(object$feature)])
    	newdata2 <- rbind(newdata,feature) 
    	newdata2 <- apply(newdata2,2,scale)
 	
    	pred <- cbind(rep(1,nrow(newdata2)),newdata2)%*%t(object$coefficient)
    	}
    
    k <- object$classes
    c <- -(1+sqrt(k))/((k-1)^(3/2))
    d <- sqrt(k/(k-1));
    vertex <- matrix(rep(0,k*(k-1)),nrow=k-1,ncol=k)
    vertex[,1] <- (k-1)^(-1/2)
    for (kk in 2:k){
      vertex[,kk] <- c
      vertex[kk-1,kk] <- c+d;
    }
    #y2 <- diag(k)%*%vertex;
    
    norm <- function(x)sqrt(sum(x^2))
    ddd <- numeric();
    for (kk in 1:k){
      ddd <- cbind(ddd,apply(scale(pred,vertex[,kk],T),1,norm))
    }
    class.pred <- apply(ddd,1,which.min)
    class.pred.out <- class.pred[1:nrow(newdata)]
    return(class.pred.out)

  }
}
