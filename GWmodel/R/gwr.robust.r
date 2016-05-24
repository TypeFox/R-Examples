gwr.robust<-function(formula, data, regression.points, bw,filtered=FALSE, kernel = "bisquare", adaptive = FALSE, p = 2, theta = 0, longlat = F, dMat, F123.test = F,maxiter=20,cut.filter= 3,cut1=2,cut2=3,delta=1.0e-5) 
    {    	   
    	W.vect<-NULL 	
    	res1<- gwr.basic(formula, data, regression.points, bw, kernel, adaptive, p, theta, longlat, dMat, F123.test, T,W.vect)
    	if(filtered==TRUE){
    		W.vect<-as.numeric(abs(res1$SDF$Stud_residual)<cut.filter)
    		
    	  	res1<-gwr.basic(formula, data, regression.points, bw, kernel, adaptive, p, theta, longlat, dMat, F123.test, T,W.vect)
    	}
    	
    	else{ # Automatic approach
    		
    		  filt <- function(x)
     		{
     			result <- rep(1,length(x))
    			xg2 <- x > cut1
      			xg3 <- x > cut2
      			span <- cut2 - cut1
      			result[xg2] <- (1 - ((x[xg2]-cut1)/span)^2)^2
      			result[xg3] <- rep(0,sum(xg3))
      			result
      		}
    		
    		iter<-0
    		diffmse<-1
    		err<-res1$SDF$residual  
    		mse<- sum(err*err)/length(err)
    		W.vect <- filt(abs(err/sqrt(mse)))
    		W.vect[is.na(W.vect)]<-0
    
    		while(diffmse>delta & iter<maxiter) {
   		 	old.mse<-mse
			res1<- gwr.basic(formula, data, regression.points, bw, kernel, adaptive, p, theta, longlat, dMat, F123.test, T,W.vect)
    	
    			err<-res1$SDF$residual 
    			mse<- sum(err*err)/length(err)
    			W.vect <- filt(abs(err/sqrt(mse)))
    			W.vect[is.na(W.vect)]<-0
    			diffmse<-abs(old.mse-mse)/mse
    			# print(round(c(iter,diffmse),3))
    			iter<-iter+1
    		}
    	}
    	
    	res1
 }
