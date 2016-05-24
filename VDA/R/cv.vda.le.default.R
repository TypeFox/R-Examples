cv.vda.le.default <-
function(x,y,kfold,lam.vec.1,lam.vec.2)
{
  if (length(y)!=nrow(x))
  stop("Dimention doesn't match! 
         Rows of feature matrix X must be the number of cases")
         
  if ((missing(lam.vec.1))||(missing(lam.vec.2)))
    stop("lambda vectors cannot be missing")
  
  if (missing(kfold))
    stop("kfold cannot be missing")
      
  y<-as.data.frame(y)
  xy<-cbind(x,y)
  xyran<-xy[sample(nrow(xy),nrow(xy)),]
  x<-xyran[,1:(ncol(xyran)-1)]
  y<-xyran[,ncol(xyran)]
  n<-length(lam.vec.1)
  m<-length(lam.vec.2)
  lam.error<-matrix(rep(0),nrow=n,ncol=m)
  lam.opt<-numeric()
  error.cv <- matrix(0,nrow=n,ncol=m)
  for(i in 1:n)
  {
    for(j in 1:m)
    {
    error.fold<-matrix(0,nrow=1,ncol=kfold)
    	for (fold in 1:kfold){
    		
   			 if (fold<kfold){
   			   n.test <- floor(length(y)/kfold)
   			   ind.data.test <- ((fold-1)*n.test+1):(fold*n.test)
  			    }
  			 else{
      			n.test <- length(y)-floor(length(y)/kfold)*(kfold-1)
     			 ind.data.test <- (length(y)-n.test+1):(length(y))
  				  }
  			n.train <- length(y)-n.test
   			ind.data.train <- which((1:length(y)) %in% ind.data.test==FALSE)
    		data.test.x <- x[ind.data.test,]
    		data.test.y <- y[ind.data.test]
    		data.train.x <- x[ind.data.train,]
    		data.train.y <- y[ind.data.train]
    
        	vda.out.train <- vda.le(data.train.x,data.train.y,lam.vec.1[i],lam.vec.2[j])
        	class.pred.test <- predict(vda.out.train,data.test.x)
        	error.fold[fold] <- length(which(as.double(class.pred.test)!=data.test.y))/length(data.test.y)
 		 							#}   
 				 }
 	   error.cv[i,j]<-sum(error.fold)/kfold
	   }
   }
  lam.min<-which(error.cv == min(error.cv), arr.ind = TRUE)
  colnames(error.cv)<-c(lam.vec.2)
  rownames(error.cv)<-c(lam.vec.1)
  
  # output
  lam.opt<-c()
  lam1.min<-lam.vec.1[lam.min[nrow(lam.min),1]]
  lam2.min<-lam.vec.2[lam.min[nrow(lam.min),2]]
  lam.opt<-c(lam1.min, lam2.min)
  est.opt<-(vda.le(x,y,lam.opt[[1]],lam.opt[[2]]))$coefficient
  error.cv<-error.cv
  
  out <- list(kfold = kfold, lam.vec.1 = lam.vec.1, lam.vec.2=lam.vec.2, error.cv = error.cv, lam.opt = lam.opt)
  class(out) <- "cv.vda.le"
  out
}
