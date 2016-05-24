bichon_optim <- function(x, model, T, method.param=NULL){
	
	y <- t(x)
  if((nrow(y)==1) && (ncol(y)==model@d)){
    z <- y #only one point
  }else{
    #many points
    if(ncol(x) == model@d) z <- x
    if(ncol(x) != model@d) z <- y
  }
	
	krig <- predict_nobias_km(object=model,newdata=as.data.frame(z),type="UK",se.compute=TRUE)
	if(is.null(method.param)) method.param <- 1
	mk    <- krig$mean; sk    <- krig$sd; alpha <- method.param
	t <- (mk-T)/sk; tplus <- t + alpha; tminus <- t - alpha

	G <- alpha*(pnorm(tplus)-pnorm(tminus)) - t*(2*pnorm(t) - pnorm(tplus) - pnorm(tminus)) - (2*dnorm(t) - dnorm(tplus) - dnorm(tminus))
			
	crit <- G*sk

	return(crit)
}

