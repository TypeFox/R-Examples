riley <- function(X, type="effect.size", optimization = "Nelder-Mead", control = list(), ...) UseMethod("riley")

riley.default <- function(X, type="effect.size", optimization = "Nelder-Mead", control = list(), ...)
{
	est <- NA    
	if (type=="test.accuracy") { est <- rileyDA(X, optimization = optimization, control=control, ...) }
	else if (type=="effect.size") { est <- rileyES(X, optimization = optimization, control=control, ...)}
	else stop(paste("Unknown type '",type,"' of meta-analysis",sep="")) 

	class(est) <- "riley"
	est
}

# effect sizes data meta-analysis
rileyES <- function(X = NULL, Y1, Y2, vars1, vars2, optimization = "Nelder-Mead", control = list(),...)
{
	if(!is.null(X)){
		X <- as.data.frame(X)
		origdata <- X
		Y1 <- X$Y1
		Y2 <- X$Y2
		vars1 <- X$vars1
		vars2 <- X$vars2
	} else {
		origdata <- cbind(Y1,vars1,Y2,vars2)
		colnames(origdata) = c("Y1","vars1","Y2","vars2")
	}
  
	numstudies = length(Y1)
	nobs <- length(which(!is.na(Y1)))+length(which(!is.na(Y2)))
	
	if(nobs != numstudies*2){warning("There are missing observations in the data!")}
	
	df <- 5 #There are 5 parameters to estimate
	if(numstudies-df < 0){warning("There are very few primary studies!")}
	
	vars = cbind(vars1, vars2)
	Y = array(NA,dim=c((length(Y1)*2),1))
	for (i in 1:length(Y1))
	{
		Y[((i-1)*2+1)] = Y1[i]
		Y[((i-1)*2+2)] = Y2[i]
	}
	
	#Calculate starting values for optim
	pars.start = c(0,0,0,0,0)
	if (numstudies >= 2) {
		sumlY1 <- uvmeta(r=Y1, vars=vars1, method="MOM")
		sumlY2 <- uvmeta(r=Y2, vars=vars2, method="MOM")
		pars.start = c(sumlY1$results["mu","Estimate"],sumlY2$results["mu","Estimate"],sqrt(sumlY1$results["mu","Var"]),sqrt(sumlY2$results["mu","Var"]),0)
	}
	
	negfullloglik <- function(pars,Y,vars)
	{
		beta1 = pars[1]
		beta2 = pars[2]
		psisq1 = pars[3]**2 #ensure variance is positive
		psisq2 = pars[4]**2 #ensure variance is positive
		rho = inv.logit(pars[5])*2-1 #ensure correlation is in [-1,1], and values in that interval move symmetric from -1 to 0 and from 1 to 0
		k = 2 #2 endpoints
		n = dim(Y)[1]/2
	
		#Beta vector
		Beta = rbind(beta1,beta2)
	
		#Design matrix
		X = array(NA,dim=c(n*2,2))
		X[,1] = rep(c(1,0),n)
		X[,2] = rep(c(0,1),n)
	
		#Create Phi matrix
		Phi = array(0,dim=c((n*2),(n*2)))
		for (i in 1:n) {
			Phi[((i-1)*2+1),((i-1)*2+1)] = vars[i,1]+psisq1
			Phi[((i-1)*2+2),((i-1)*2+2)] = vars[i,2]+psisq2
			Phi[((i-1)*2+1),((i-1)*2+2)] = rho*sqrt((vars[i,1]+psisq1)*(vars[i,2]+psisq2))
			Phi[((i-1)*2+2),((i-1)*2+1)] = rho*sqrt((vars[i,1]+psisq1)*(vars[i,2]+psisq2))
		}
		
		#Minimize the negative of the restricted log-lkh
		0.5*((n-k)*log(2*pi)-log(det(t(X)%*%X))+log(det(Phi))+log(det(t(X)%*%solve(Phi)%*%X))+(t(Y-X%*%Beta)%*%solve(Phi)%*%(Y-X%*%Beta)))
	}
	
	fit = optim(pars.start,negfullloglik,Y=Y,vars=vars,method=optimization,hessian=T,control=control)
	
	if(fit$convergence != 0) { 
		if(fit$convergence == 1) warning ("Iteration limit had been reached.")
		else if (fit$convergence == 10) warning("Degeneracy of the Nelder-Mead simplex.")
		else if (fit$convergence == 51 | fit$convergence == 52) warning(fit$message)
    else warning("Unspecified convergence error in optim.")
	}
	
	beta1 = fit$par[1]
	beta2 = fit$par[2]
	psi1 = abs(fit$par[3])
	psi2 = abs(fit$par[4])
	rhoT = fit$par[5]
	coefficients = c(beta1,beta2,psi1,psi2,rhoT)
	names(coefficients) = c("beta1","beta2","psi1","psi2","rhoT")
	
	hessian = fit$hessian
	colnames(hessian) = c("beta1","beta2","psi1","psi2","rhoT")
	rownames(hessian) = c("beta1","beta2","psi1","psi2","rhoT")
	
	if (length(which(eigen(fit$hessian,symmetric=TRUE)$values<0))>0) warning("The Hessian contains negative eigenvalues!")
	
	iterations <- fit$iterations
	logLik <- -fit$value
	
	output <- list(coefficients = coefficients, hessian = hessian, df = df, numstudies = numstudies, nobs = nobs, logLik = logLik,
			   iterations = (iterations+1), call = match.call(), data = origdata, type="effect.size")  
	return(output)
}

# Diagnostic test accuracy data meta-analysis
rileyDA <-
  function(X = NULL, TP, FN, FP, TN, correction = 0.5, 
           correction.control = "all", optimization = "Nelder-Mead", control = list(), ...)
  {
      if(!is.null(X)){
        X <- as.data.frame(X)
        origdata <- newdata <- X
      } else {
        origdata <- newdata <- as.data.frame(cbind(TP,FN,FP,TN))
        colnames(origdata) <- c("TP","FN","FP","TN")
      }
      
	  ## The following corrections are copied from the "mada" package to facilitate comparison of results
      ## apply continuity correction to _all_ studies if one contains zero
      if(correction.control == "all"){if(any(origdata == 0)){newdata$TP <- origdata$TP + correction;
                                                             newdata$FN <- origdata$FN + correction;
                                                             newdata$FP <- origdata$FP + correction;
                                                             newdata$TN <- origdata$TN + correction}}
      if(correction.control == "single"){
        correction = ((((origdata$TP == 0)|(origdata$FN == 0))|(origdata$FP == 0))| (origdata$TN == 0))*correction
        newdata$TP <- correction + origdata$TP
        newdata$FN <- correction + origdata$FN
        newdata$FP <- correction + origdata$FP
        newdata$TN <- correction + origdata$TN
      }
      
      
      #Calculate sensitivities and specificities (original scale)
      number.of.pos <- newdata$TP + newdata$FN
      number.of.neg <- newdata$FP + newdata$TN
      sens <-newdata$TP/number.of.pos
      fpr <- newdata$FP/number.of.neg
      var.sens = sens*(1-sens)/number.of.pos
      var.fpr = fpr*(1-fpr)/number.of.neg
      
      logit.sens <- logit(sens)
      logit.fpr <- logit(fpr)
      var.logit.sens <- 1/(sens*(1-sens)*number.of.pos)
      var.logit.fpr <- 1/(fpr*(1-fpr)*number.of.neg)
      
	    #Apply ordinary bivariate meta-analysis on transformed data
      output = rileyES(X=NULL, Y1=logit.sens,Y2=logit.fpr,vars1=var.logit.sens,vars2=var.logit.fpr,optimization = optimization, control = control, ...)
      output$type = "test.accuracy"
      output$data = newdata
      output$correction = correction 
      output$correction.control = correction.control
      
      return(output)
  }

print.riley <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  if (length(which(eigen(x$hessian,symmetric=TRUE)$values<0))>0) cat("\nWarning: the Hessian matrix contains negative eigenvalues, parameter estimates are thus not optimally fitted!\n")
}

# Calculate prediction interval (not identical interpretation to random effects!)
predict.riley <- function(object, level = 0.95, ...)
{
  alpha <- (1-level)/2
  
  predint		<- array(NA,dim=c(2,3))
  colnames(predint) <- c("Estimate", paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  df 			<- object$df
  numstudies 		<- object$numstudies
  
  if (object$type=="test.accuracy")
  {
	rownames(predint) = c("Sens","FPR")
	if ((numstudies - df) > 0)
	{
		predint[1,] = inv.logit(c(coefficients(object)["beta1"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"])))
		predint[2,] = inv.logit(c(coefficients(object)["beta2"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"])))
	} else {
		predint[1,] = c(inv.logit(coefficients(object)["beta1"]),0,1)
		predint[2,] = c(inv.logit(coefficients(object)["beta2"]),0,1)
	}
  } else {
	rownames(predint) = c("beta1","beta2")
	predint[1,] = c(coefficients(object)["beta1"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi1"]**2)+diag(vcov(object))[1])+coefficients(object)["beta1"]))
	predint[2,] = c(coefficients(object)["beta2"],(qt(alpha,df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]),(qt((1-alpha),df=(numstudies-df))*sqrt((coefficients(object)["psi2"]**2)+diag(vcov(object))[2])+coefficients(object)["beta2"]))
  }
  predint
}


summary.riley <- function(object, level = 0.95, ...)
{
	confints <- cbind(object$coefficients, confint(object,level=level))
	colnames(confints)[1] <- "Estimate"
	
	if (object$type=="test.accuracy") {
		confints <- rbind(inv.logit(confints[1:2,]),confints)
		rownames(confints)[1:2] <-  c("Sens", "FPR") 
	} 
	
	#Transform last parameter back to rho
	confints["rhoT",] =  inv.logit(confints["rhoT",])*2-1
	rownames(confints)[which(rownames(confints)=="rhoT")] = "rho"
	
	res <- list(call=object$call, confints = confints)
	class(res) <- "summary.riley"
	res
}

print.summary.riley <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$confints)
}

vcov.riley <- function(object, ...){
  if (length(which(eigen(object$hessian,symmetric=TRUE)$values<0))>0) warning("The Hessian contains negative eigenvalues!")
 
  # It is known that 'optim' has problems.  Perhaps the simplest thing to do is to call 'optim' with each of the 'methods' in sequence, using the 'optim' found by each 'method' as the starting value for the next.  When I do this, I often skip 'SANN', because it typically takes so much more time than the other methods.  However, if there might be multiple local minima, then SANN may be the best way to find a global minimum, though you may want to call 'optim' again with another method, starting from optimal solution returned by 'SANN'. 

  Sigma = solve(object$hessian)
  Sigma
  }

logLik.riley <- function(object, ...){
	val 				<- object$logLik
	attr(val, "nobs") 	<- object$nobs
	attr(val, "df") 	<- object$df
	class(val) 			<- "logLik"
	val
}


plot.riley <- function(x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, main="",
                       ylim = c(0,1), xlim = c(0,1), pch = 1, lty = 1, lwd = 1, cex.numerics=0.45,
                       add=FALSE, ...)
{
	alpha = (1-level)/2
	
	if (x$type=="test.accuracy") {
		xlab = "1-Specificity"
		ylab = "Sensitivity"
		
		FP <- x$data$FP
		negatives <- FP + x$data$TN
		FPR <- FP/negatives
		mu = x$coefficients[c("beta2","beta1")]
		Sigma = vcov(x)[c("beta2","beta1"),c("beta2","beta1")] 
		mu.ellipse <- ellipse(Sigma, centre = mu, level = level) 
		summary1 = inv.logit(mu[1])
		summary2 = inv.logit(mu[2])
		ellipse1 = inv.logit(mu.ellipse[,1])
		ellipse2 = inv.logit(mu.ellipse[,2])
	} else {
		plotnumerics = FALSE
		xlab = "Y1"
		ylab = "Y2"
		mu = x$coefficients[c("beta1","beta2")]
		Sigma = vcov(x)[c(1,2),c(1,2)] 
		mu.ellipse <- ellipse(Sigma, centre = mu, level = level) 
		summary1 = mu[1]
		summary2 = mu[2]
		ellipse1 = mu.ellipse[,1]
		ellipse2 = mu.ellipse[,2]
	}
	
	if (!add) plot(-500,-500, type = "l", xlim = xlim, ylim = ylim, xlab=xlab,ylab=ylab,main=main, ...)
	#if (!add) NextMethod("plot")	
	polygon(ellipse1,ellipse2,lty=lty, lwd=lwd)
	if(plotsumm) points(summary1,summary2,pch=pch) # add the point estimate of the mean
	
	if(plotnumerics) {
		ci = summary(x,level=level)[2]$confints
		text(0.8,0.15,labels="Estimate",pos=2,cex=cex.numerics)
		text(0.9,0.15,labels=paste((alpha*100),"% CI",sep=""),pos=2,cex=cex.numerics)
		text(1.0,0.15,labels=paste(((1-alpha)*100),"% CI",sep=""),pos=2,cex=cex.numerics)
		text(0.5,0.10,labels= "Sensitivity",pos=4, cex=cex.numerics)
		text(0.8,0.10,labels=paste("",formatC(round( ci["Sens",1],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(0.9,0.10,labels=paste("",formatC(round( ci["Sens",2],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(1.0,0.10,labels=paste("",formatC(round( ci["Sens",3],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(0.5,0.05,labels= "Specificity",pos=4, cex=cex.numerics)
		text(0.8,0.05,labels=paste("",formatC(round( 1-ci["FPR",1],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(0.9,0.05,labels=paste("",formatC(round( 1-ci["FPR",3],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
		text(1.0,0.05,labels=paste("",formatC(round( 1-ci["FPR",2],2),2,format="f",flag="0")),pos=2, cex=cex.numerics)
	}
}
