# source("EI.R")
EI <- function (x, model, plugin=NULL, type="UK", minimization = TRUE, envir=NULL) {
  
   if (is.null(plugin)){ 
    if (minimization) {
      plugin <- min(model@y)
    } else {
      plugin <- -max(model@y)
    }
   }
	 m <- plugin

   ########################################################################################
   # Convert x in proper format(s)
   d <- length(x)
   if (d != model@d){ stop("x does not have the right size") }
   newdata.num <- as.numeric(x)
   newdata <- data.frame(t(newdata.num))
   colnames(newdata) = colnames(model@X)

   ########################################################################################
   predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE)
   kriging.mean <- predx$mean
   if(!minimization) {
    kriging.mean <- -kriging.mean
   }
   kriging.sd   <- predx$sd

   xcr <- (m - kriging.mean)/kriging.sd
    
	if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) 
	{ res <- 0
    xcr <- xcr.prob <- xcr.dens <- NULL
	} else 
  {   xcr.prob <- pnorm(xcr)
      	xcr.dens <- dnorm(xcr)	        
	   	  res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
	}
    
  if (!is.null(envir)) 
  { assign("xcr", xcr, envir=envir)
    assign("xcr.prob", xcr.prob, envir=envir)
    assign("xcr.dens", xcr.dens, envir=envir)
    assign("kriging.sd", kriging.sd, envir=envir)
		assign("c", predx$c, envir=envir)
 		assign("Tinv.c", predx$Tinv.c, envir=envir)
    
 	}
	return(res)
}
