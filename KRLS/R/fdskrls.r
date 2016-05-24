fdskrls <- 
  function(object,...){ 
    d <- ncol(object$X)
    n <- nrow(object$X)
    lengthunique    <- function(x){length(unique(x))}
    
    fdderivatives           =object$derivatives
    fdavgderivatives        =object$avgderivatives
    fdvar.var.avgderivatives=object$var.avgderivatives    
    
    # vector with positions of binary variables
    binaryindicator <-which(apply(object$X,2,lengthunique)==2)    
      if(length(binaryindicator)==0){
        # no binary vars in X; return derivs as is  
      } else {
        # compute marginal differences from min to max 
        est  <- se <- matrix(NA,nrow=1,ncol=length(binaryindicator))
        diffsstore <- matrix(NA,nrow=n,ncol=length(binaryindicator))
        for(i in 1:length(binaryindicator)){
          X1 <- X0 <- object$X
          # test data with D=Max
          X1[,binaryindicator[i]] <- max(X1[,binaryindicator[i]])
          # test data with D=Min
          X0[,binaryindicator[i]] <- min(X0[,binaryindicator[i]])
          Xall      <- rbind(X1,X0)
          # contrast vector
          h         <- matrix(rep(c(1/n,-(1/n)),each=n),ncol=1)
          # fitted values
          pout      <- predict(object,newdata=Xall,se=TRUE)
          # store FD estimates
          est[1,i] <- t(h)%*%pout$fit        
          # SE (multiply by sqrt2 to correct for using data twice )
          se[1,i] <- as.vector(sqrt(t(h)%*%pout$vcov.fit%*%h))*sqrt(2)
          # all
          diffs <- pout$fit[1:n]-pout$fit[(n+1):(2*n)]          
          diffsstore[,i] <- diffs 
        }        
        # sub in first differences
        object$derivatives[,binaryindicator] <- diffsstore
        object$avgderivatives[,binaryindicator] <- est
        object$var.avgderivatives[,binaryindicator] <- se^2
        object$binaryindicator[,binaryindicator] <- TRUE
      }  
    
   return(invisible(object))
  
}
  