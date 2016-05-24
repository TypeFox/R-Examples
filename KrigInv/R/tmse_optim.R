tmse_optim <- function(x, model, T, method.param=NULL){

	y <- t(x)
	if(ncol(y) == model@d) z <- y
	if(ncol(y) != model@d) z <- x
	
	krig <- predict_nobias_km(object=model,newdata=as.data.frame(z),type="UK",se.compute=TRUE)
	mk    <- krig$mean
	sk    <- krig$sd
	
	if(is.null(method.param)) method.param <- 0
	epsilon <- method.param
	W <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2+epsilon^2))^2)
	tmse <- W*sk^2
	return(tmse)
}
