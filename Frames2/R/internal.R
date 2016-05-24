Lag1 <- function (u, ds, mu)
{
	L <- -1 / max(u - mu)
	R <- -1 / min(u - mu)
	dif <- 1
	tol <- 1e-04
	while (dif > tol){
		M <- (L + R) / 2
		glam <- sum((ds * (u - mu)) / (1 + M * (u - mu)))
		if (glam > 0) L <- M
		if (glam < 0) R <- M
		dif <- abs (glam)
	}
	return (M)
}

Lag2 <- function (u, ds, mu)
{
	n <- length (ds)
	u <- u - rep (1,n) %*% t(mu)
	M <- 0 * mu
	dif <- 1
	tol <- 1e-04
	while (dif > tol){
		D1 <- 0 * mu
		DD <- D1 %*% t(D1)
		for (i in 1:n){
			aa <- as.numeric (1 + t(M) %*% u[i,])
			D1 <- D1 + ds[i] * u[i,] / aa
			DD <- DD - ds[i] * (u[i,] %*% t(u[i,]))/aa^2
		}
		D2 <- solve (DD, D1, tol = 1e-12)
		dif <- max (abs(D2))
		rule <- 1
		while (rule > 0){
			rule <- 0
			if (min (1 + t(M - D2) %*% t(u)) <= 0) rule <- rule + 1
			if (rule > 0) D2 <- D2 / 2
		}
		M <- M - D2	
	}
	return (M)
}

PELConfInt = function (a, ys, ds, YEL, nss, ps) 
{
	res <- rep (0, 2)
	#--------------------------------------
	tol <- 1e-08
	cut <- qchisq(a,1)
	#--------------
	t1 <- YEL
	t2 <- max(ys)
	dif <- t2 - t1
	while (dif > tol) {
		tau <- (t1 + t2) / 2
		M <- Lag1 (ys, ds, tau)
		elratio <- -2 * nss * (sum(ds * (log(ds) - log (1 + M * (ys - tau)) - log(ps))))
		if (elratio > cut) t2 <- tau
		if (elratio <= cut) t1 <- tau
		dif <- t2 - t1
	}
	res[1] <- (t1 + t2) / 2
	#-------------
	t1 <- YEL
	t2 <- min(ys)
	dif <- t1 - t2
	while (dif > tol) {
		tau <- (t1 + t2) / 2
		M <- Lag1(ys, ds, tau)
		elratio <- -2 * nss * (sum(ds * (log(ds) - log (1 + M * (ys - tau)) - log(ps))))
		if(elratio > cut) t2 <- tau
		if(elratio <= cut) t1 <- tau
		dif <- t1 - t2
	}
	res[2] <- (t1 + t2) / 2
	return(res)
}

`print.EstimatorDF`=function(x, ...){
  if (is.null(attr(x, "attributesDF")[1])){
  	cat("\nEstimation:\n")
  	print(x$Est)
  }
  else {
	cat("\nEstimation and ",attr(x, "attributesDF")[1]*100,"% Confidence Intervals:\n")
  	print(x$ConfInt)
  }
}

`summary.EstimatorDF`=function(object, ...){
  cat("\nCall:\n")
  print(object$Call)
  cat("\nEstimation:\n")
  print(object$Est)
  if (!is.null(object$VarEst)){
	cat("\nVariance Estimation:\n")
  	print(object$VarEst)
  }
  if (!is.null(object$TotDomEst)){
  	cat("\nTotal Domain Estimations:\n")
  	print(object$TotDomEst)
  }
  if (!is.null(object$MeanDomEst)){
  	cat("\nMean Domain Estimations:\n")
  	print(object$MeanDomEst)
  }
  if(!is.null(object$Param)){
  	cat("\nParameters:\n")
  	print(object$Param)
  }
  if (!is.null(attr(object, "attributesDF")[1])){
 	cat("\n",attr(object, "attributesDF")[1]*100,"% Confidence Intervals:\n")
  	print(object$ConfInt)
  }
}

`print.EstimatorMDF`=function(x, ...){
  if (is.null(attr(x, "attributesMDF")[1])){
  	cat("\nEstimation:\n")
  	print(x$Est)
  }
  else {
	cat("\nEstimation and ",attr(x, "attributesMDF")[1]*100,"% Confidence Intervals:\n")
  	print(x$ConfInt)
  }
}
