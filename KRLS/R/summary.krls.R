summary.krls <-
function(object, probs=c(.25,.5,.75),...)
      {
            
        if( class(object)!= "krls" ){
        warning("Object not of class 'krls'")
        UseMethod("summary")
        return(invisible(NULL))
        }
        
        cat("* *********************** *\n")
        cat("Model Summary:\n\n")
        cat("R2:",object$R2,"\n\n")
        
        d <- ncol(object$X)
        n <- nrow(object$X)
        
        coefficients <- matrix(NA,d,0)
        rownames(coefficients) <- colnames(object$X) 
                
       
        if(is.null(object$derivatives)){
          cat("\n")
          cat("recompute krls object with krls(...,derivative = TRUE) to get summary of marginal effects\n")
          return(invisible(NULL))
        } 
        
        # average marginal effects  
        est     <- t(object$avgderivatives)
        se     <- sqrt(t(object$var.avgderivatives))
        tval   <- est/se
        avgcoefficients <- cbind(est, se, tval, 2 * pt(abs(tval),n-d, lower.tail = FALSE))
        colnames(avgcoefficients) <- c("Est", "Std. Error", "t value", "Pr(>|t|)")
        
       # add stars for binary    
        if(sum(object$binaryindicator)>0){         
          rownames(avgcoefficients)[object$binaryindicator] <- paste(rownames(avgcoefficients)[object$binaryindicator],"*",sep="")
        }
        
        cat("Average Marginal Effects:\n")
        print(avgcoefficients,...)
        if(sum(object$binaryindicator)>0){
        cat("\n(*) average dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
        }
        
        # quantiles of derivatives
        qderiv <- apply(object$derivatives,2,quantile,probs=probs)
        if(sum(object$binaryindicator)>0){         
          colnames(qderiv)[object$binaryindicator] <- paste(colnames(qderiv)[object$binaryindicator],"*",sep="")
        }
        qderiv <- t(qderiv)
        
        cat("\n")
        cat("Quartiles of Marginal Effects:\n")
        print(qderiv,...)
        
        if(sum(object$binaryindicator)>0){
         cat("\n(*) quantiles of dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
        }
                     
      ans <- list(
                 coefficients=avgcoefficients,
                 qcoefficients=qderiv)
      class(ans) <- "summary.krls"  
      return(invisible(ans))
}

