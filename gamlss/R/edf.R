#======================================================================================
edf <- function(obj, what = c("mu", "sigma", "nu", "tau"), parameter= NULL,  print=TRUE, ...)
 {
what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
 x <- obj[[paste(what,"coefSmo",sep=".")]]
 s <- obj[[paste(what,"s",sep=".")]]
 P <- length(x)
 edf  <- rep(0, length=P)
 for (i in 1:P)
 edf[i] <- x[[i]]$edf
  names(edf) <- colnames(s)
 if (print) cat("Effective df for", what, "model", "\n")
   edf
 } 
#----------------------------------------------------------------------------------------
edfAll <-function(obj,  ...)       
 {
      out <- list()
     if ("mu" %in% obj$par) #
         out$mu <- edf(obj, "mu" , print=FALSE)
     if ("sigma" %in% obj$par)  
      out$sigma <- edf(obj, "sigma", print=FALSE )
     if ("nu" %in% obj$par)  
         out$nu <- edf(obj, "nu", print=FALSE )
     if ("tau" %in% obj$par)  
        out$tau <-  edf(obj, "tau", print=FALSE)
       return(out)
 }
#----------------------------------------------------------------------------------------    
  
