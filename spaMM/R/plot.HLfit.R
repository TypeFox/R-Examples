`plot.HLfit` <- function(x, which=c("mean","ranef"),
                         titles = list(
                           meanmodel=list(outer="Mean model",devres="Deviance residuals", absdevres="|Deviance residuals|",
                                          resq="Residual quantiles", devreshist="Deviance residuals"),
                           ranef=list(outer="Random effects and leverages",qq="Random effects Q-Q plot",
                                      levphi=expression(paste("Leverages for ",phi)), levlambda=expression(paste("Leverages for ",lambda)))
                           ),
                         control = list() , ...) {
	residuals <- x$std_dev_res 
 	fitted.values <- x$fv
  ## possible modif of 'which':
  if (is.null(residuals)) {## possible if disp pars in ranFix
    which <- setdiff(which,c("mean")) ## all 'mean' diagnostic plots involve these residuals
  }
  if ("predict" %in% which) {
    plot(x$y,predict(x),xlab="Response",ylab="Predicted response",...)
    abline(0,1)
  }
	if ("ranef" %in% which) {
	  if (x$family$family %in% c("poisson","binomial","COMPoisson")) {
	    lev_phi <- NULL  ## currently (10/2013) <HLfit>$lev_phi non-null even for poisson, binomial
	  } else lev_phi <- x$lev_phi
	  lev_lambda <- x$lev_lambda
	  ranef <- x$ranef ## u
	  nranplots <- length(c(which(length(lev_lambda)>0),which(length(lev_phi)>0),which(length(ranef)>0)))
    if (nranplots==0L) which <- which[which!="ranef"]
	} 
  pch <- control$pch  
  if (is.null(pch)) pch <- "+"   
  pcol <- control$pcol  
  if (is.null(pcol)) pcol <- "blue"  
  lcol <- control$lcol  
  if (is.null(lcol)) lcol <- "red" 
  if (interactive()) par(ask=TRUE) 
  for (i in seq_len(length(which))) {
    typ <- which[i]
    if (i > 1) dev.new() ## =: meaning of having different elements in 'which'  
		if (typ =="mean") { ## diagnostic plots for mean model => 4 subplots
		  par(mfrow = c(2, 2), oma = c( 0, 0, 2, 0 ), pty = "s", ...) ## 4 subplots for mean !!
      #
		  loess.fit <- loess.smooth(fitted.values, residuals)
			plot(fitted.values, residuals, xlab = "Fitted Values", 
				 ylab = titles$meanmodel$devres, pch = pch, col = pcol, bty = "n", main = titles$meanmodel$devres)
			lines(loess.fit$x, loess.fit$y, col = lcol)
			#
      loess.fit <- loess.smooth(fitted.values, abs(residuals))
			plot(fitted.values, abs(residuals), xlab = "Fitted Values", 
				 ylab = titles$meanmodel$absdevres, pch = pch, col = pcol, bty = "n", main = titles$meanmodel$absdevres)
			lines(loess.fit$x, loess.fit$y, col = lcol)
			#
      qqnorm(residuals, col = pcol, pch = pch, bty = "n", 
				   xlab = "Normal quantiles", ylab = titles$meanmodel$resq, main = titles$meanmodel$resd)
			qqline(residuals, col = lcol)
			#
      hist(residuals, density = 15, xlab = titles$meanmodel$devreshist, main = "", col = pcol)
            title(titles$meanmodel$outer,outer=TRUE)
		}
		if (typ == "ranef") {
		  if (nranplots<3) {
		    par(mfrow = c(1, nranplots), oma = c( 0, 0, 2, 0 ), pty = "s", ...)
		  } else { par(mfrow = c(2, 2), oma = c( 0, 0, 2, 0 ), pty = "s", ...) } ## 
		  if (length(ranef)>0) {
              std_ranef <- ranef/sqrt(1-lev_lambda)
     		  qqnorm(ranef, col = pcol, pch = pch, bty = "n", 
				   xlab = "Normal quantiles", ylab = expression(paste("Standardized ",italic(u)," quantiles")), main = titles$ranef$qq)
	    	  qqline(ranef, col = lcol)
            }
			if (!is.null(lev_phi)) plot(lev_phi, ylab = "", main = titles$ranef$levphi , pch = pch, col = pcol, bty = "n")
			if (length(lev_lambda)>0) plot(lev_lambda, ylab = "", main= titles$ranef$levlambda , pch = pch, col = pcol, bty = "n")
		  title(titles$ranef$outer,outer=TRUE)
		} 
  }
	par(ask=FALSE) 
	invisible(x)
}


