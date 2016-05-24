robpredict <-
function(fit, areameans=NULL, k=NULL, reps=NULL){
   if (!inherits(fit, "fitsaemodel")) stop("fit must be of class 'fitsaemodel'")
   # get the modelk
   modelk <- attr(fit, "method")$tuning$k
   # get the decomposition
   dec <- attr(fit, "dec")
   # for ml, modelk = 20000
   if (is.null(modelk)) modelk <- 20000
   # prediction k (which is not necessarily equal to modelk) 
   if (is.null(k)){
      # here: k == modelk
      k <- modelk 
   }else{
      if (k <= 0 ) stop("Robustness tuning constant k must be > 0!\n") 
   }
   kappa <- .computekappa(k)
   # initialize
   mspe <- NULL
   # from the model definitions; used in order to compute the random effects
   model <- attr(fit, "saemodel")
   areaNames <- attr(model, "areaNames")
   x <- model$X
   y <- model$y
   n <- model$n
   p <- model$p
   g <- model$g
   nsize <- model$nsize
   # from the fitted model
   beta <- fit$beta
   v <- fit$theta[1]
   d <- fit$theta[2] / v
   # tau
   tau <- c(beta, v, d)
   # preparations for fortran-call
   predre <- rep(0, g)
   predfe <- rep(0, g)
   tmp <- .Fortran("drsaehubpredict", n=as.integer(n), p=as.integer(p), g=as.integer(g), nsize=as.integer(nsize), k=as.double(k), kappa=as.double(kappa), d=as.double(d), v=as.double(v), beta=as.matrix(beta), yvec=as.matrix(y), xmat=as.matrix(x), predfe=as.matrix(predfe), predre=as.matrix(predre), dec=as.integer(dec))
   # retrieve the area-level random effects; it is used whether new data is present or not
   raneff <- tmp$predre
   # branch: old vs new data
   if (is.null(areameans)){
      fixeff <- as.matrix(tmp$predfe)
   }else{
      # check whether the new data are proper
      if(!is.matrix(areameans)){
	 areameans <- as.matrix(areameans)
      }
      # check the dimensions
      if (dim(areameans)[1] != g) stop("'areameans' is not of conformable size! \n")
      if (dim(areameans)[2] != p) stop("'areameans' is not of conformable size! \n")
      # compute the fixed-effect spredictions (at the area level)
      fixeff <- areameans %*% beta
      # compute mspe, given the fitted model (this option is only valid if areameans != NULL
      if (!is.null(reps)){
	 mspe <- .mspe(fit, abs(round(reps)), areameans, fixeff)
      }
   }
   means <- raneff + fixeff
   rownames(fixeff) <- areaNames
   rownames(raneff) <- areaNames
   rownames(means) <- areaNames
   # compute the residuals of the model (i.e. e_ij = y_ij - X_ij*beta - u_i)
   vn <- numeric(n)
   getres <- .Fortran("drsaeresid", n=as.integer(n), p=as.integer(p), g=as.integer(g), nsize=as.integer(nsize), k=as.double(modelk), tau=as.matrix(tau), u=as.matrix(raneff), xmat=as.matrix(x), yvec=as.matrix(y), res=as.matrix(vn), stdres=as.matrix(vn), wgt=as.matrix(vn), dec=as.integer(dec))
   #
   result <- list(fixeff=fixeff, raneff=raneff, means=means, res=getres$res, stdres=getres$stdres, wgt=getres$wgt, mspe=mspe)
   attr(result, "robustness") <- k
   attr(result, "fit") <- fit
   attr(result, "mspe") <- reps
   class(result) <- "meanssaemodel"
   return(result)
}

