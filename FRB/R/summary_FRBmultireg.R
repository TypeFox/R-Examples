summary.FRBmultireg <- function(object, confmethod = c("BCA","basic","both"), 
digits=3, print.CI=FALSE, sep="", ...) {

confmethod <- match.arg(confmethod)
estBeta <- object$coefficients
q <- ncol(estBeta)
p <- nrow(estBeta)
responses <- colnames(as.data.frame(estBeta))
covariates <- rownames(as.data.frame(estBeta))

res=object[c("call", "coefficients","residuals","scale","Sigma",
            "weights", "df", "method", "control")]
res$responses=responses
res$covariates=covariates
res$digits=digits

if (is.null(object$bootest))  res$Beta=estBeta
else{
  leftbr <- rep("(", p)
  rightbr <- rep(")", p)
  Betawstd <- as.data.frame(list(estBeta[,1], leftbr, object$SE[,1], rightbr))
  if (q>1) {
    for (i in 2:q) {
        Betawstd <- as.data.frame(list(Betawstd, estBeta[,i], leftbr, object$SE[,i], rightbr))
    }
  }
  colnames(Betawstd)[1+(0:(q-1))*4] <- responses
  colnames(Betawstd)[2+(0:(q-1))*4] <- " "
  colnames(Betawstd)[3+(0:(q-1))*4] <- " "
  colnames(Betawstd)[4+(0:(q-1))*4] <- " "
  
  separator <- rep(sep, p)
  sigcodes.bca <- matrix(rep("",p*q), ncol=q)
  sigcodes.bca[object$p.bca<0.1] <- "."  
  sigcodes.bca[object$p.bca<0.05] <- "*"  
  sigcodes.bca[object$p.bca<0.01] <- "**"  
  sigcodes.bca[object$p.bca<0.001] <- "***"  
  sigcodes.basic <- matrix(rep("",p*q), ncol=q)
  sigcodes.basic[object$p.basic<0.1] <- "."  
  sigcodes.basic[object$p.basic<0.05] <- "*"  
  sigcodes.basic[object$p.basic<0.01] <- "**"  
  sigcodes.basic[object$p.basic<0.001] <- "***"  
  if ((confmethod == "BCA") | (confmethod == "both")) {
    llist.bca <- list()
    for (i in 1:q) {
  
#      limits <- cbind(object$CI.bca.lower[,i], object$CI.bca.upper[,i])
#      rownames(limits) <- rownames(estBeta)
#      colnames(limits) <- c("lower", "upper")
      if (print.CI) {BetaTable <- as.data.frame(list(estBeta[,i,drop=FALSE], separator, object$SE[,i], separator, object$CI.bca.lower[,i], object$CI.bca.upper[,i], separator, object$p.bca[,i], sigcodes.bca[,i]))
      colnames(BetaTable) <- c("Estimate", " ", "Std.Error", " ", "lower", "upper", " ", "p-value", " ")}
      else {BetaTable <- as.data.frame(list(estBeta[,i,drop=FALSE], separator, object$SE[,i], separator, object$p.bca[,i], sigcodes.bca[,i]))
      colnames(BetaTable) <- c("Estimate", " ", "Std.Error", " ", "p-value", " ")}
      llist.bca[[responses[i]]] <- BetaTable
    }
  }
  if ((confmethod == "basic") | (confmethod == "both")) {
    llist.basic <- list()
    for (i in 1:q) {
  
      if (print.CI) {BetaTable <- as.data.frame(list(estBeta[,i,drop=FALSE], separator, object$SE[,i], separator, object$CI.basic.lower[,i], object$CI.basic.upper[,i], separator, object$p.basic[,i], sigcodes.basic[,i]))
      colnames(BetaTable) <- c("Estimate", " ", "Std.Error", " ", "lower", "upper", " ", "p-value", " ")}
	else {BetaTable <- as.data.frame(list(estBeta[,i,drop=FALSE], separator, object$SE[,i], separator, object$p.basic[,i], sigcodes.basic[,i]))
      colnames(BetaTable) <- c("Estimate", " ", "Std.Error", " ", "p-value", " ")} 
      llist.basic[[responses[i]]] <- BetaTable
    }
  }
  res$Betawstd=Betawstd
  res$conf=object$conf
  res$print.CI=print.CI

  if (confmethod == "BCA") res$table.bca=llist.bca
  else if (confmethod == "basic") res$table.basic=llist.basic
  else
    {res$table.bca=llist.bca 
    res$table.basic=llist.basic}
}

class(res) <- "summary.FRBmultireg"

res

}

###############################################################################

print.summary.FRBmultireg <- function(x, ...) {

if (x$method$est=="MM")
cat(paste("\nMultivariate regression based on MM-estimates (bdp = ", x$method$bdp, ", eff = ", x$method$eff, ")", sep=""), "\n")
else
cat(paste("\nMultivariate regression based on ", x$method$est, "-estimates (breakdown point = ", x$method$bdp, ")", sep=""), "\n")
cat("\n")
 if (!is.null(x$call)) 
  cat("Call:\n", deparse(x$call), "\n", sep = "",collapse = "\n") 
  print.CI=x$print.CI
  resid <- x$residuals
  digits=x$digits	
  if (x$df > 5) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        if (NROW(resid) > 1) 
            rq <- structure(apply(resid, 2, quantile), dimnames = list(nam, 
                dimnames(resid)[[2]]))
        else rq <- structure(apply(resid, 1, quantile), dimnames = list(nam, 
                dimnames(resid)[[1]]))
      rq=t(rq)	
   }
   else {
       if (NROW(resid) > 1) rq=t(resid)
       else rq=resid
   }  

#if (!is.null(x$Betawstd)) {
#  cat("Coefficient estimates (with bootstrap standard errors):\n")
#  print(x$Betawstd, digits=x$digits)
#} 
#else {
#  cat("Coefficient estimates):\n")
#  print(x$coefficients, digits=x$digits)
#} 
#for (i in 1:length(x$table.bca))
for (i in 1:NROW(rq))
{ 
      if (NROW(rq)>1) cat("\nResponse ", x$responses[i], ":\n\n", sep="")
      cat(if (!is.null(x$weights) && diff(range(x$weights))) 
        "Residuals:\n", sep = "")
      print(rq[i,], digits = digits, ...)
	cat("\n")	
      cat("Coefficients:\n")
	if (!is.null(x$table.bca)) 
	{
	print(x$table.bca[[x$responses[i]]], digits = x$digits)
      cat("---\n") 
      cat("Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n\n")
    	if (print.CI) cat("Confidence limits (",x$conf*100, "%) and p-values based on BCA method!\n\n", sep="")
      else cat("p-values based on BCA method!\n\n", sep="")
	}

	if (!is.null(x$table.basic)) 
	{
      print(x$table.basic[[x$responses[i]]], digits = x$digits)
      cat("---\n") 
      cat("Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n\n")
      if (print.CI) cat("Confidence limits (",x$conf*100, "%) and p-values based on Basic Bootstrap!\n\n", sep="")  
	else cat("p-values based on Basic Bootstrap!\n\n", sep="")  
  	}
	if (is.null(x$table.bca) && is.null(x$table.basic))
	{
#  	cat("Coefficient estimates:\n")
      estim=x$coefficients[,i,drop=FALSE]
      colnames(estim)="Estimate"
	print(estim, digits=x$digits)
      cat("\n\n")
	}
} 

cat(paste("Robust residual scale:",format(signif(x$scale,digits)),"\n"))

cat("\nError covariance matrix estimate:\n")
print(as.data.frame(x$Sigma), digits=x$digits)
}

###############################################################################

print.FRBmultireg <- function(x, digits=3, ...) {


estBeta <- x$coefficients
q <- ncol(estBeta)
p <- nrow(estBeta)

 if (x$method$est=="MM")
  cat(paste("\nMultivariate regression based on multivariate MM-estimates (bdp = ", x$method$bdp, ", eff = ", x$method$eff, ")", sep=""), "\n")
  else
  cat(paste("\nMultivariate regression based on multivariate ", x$method$est, "-estimates (breakdown point = ", x$method$bdp, ")", sep=""), "\n")

 if (!is.null(x$call)) 
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\n") 
  cat("Coefficients:\n")
  if (q==1) print(t(estBeta)[1,,drop=FALSE], digits=digits)
  else print(estBeta, digits=digits)

}

###############################################################################

predict.FRBmultireg <- function (object, newdata, ...) 
{
    class(object) <- c(class(object), "mlm")
    object$qr <- qr(sqrt(object$weights) * object$X)
    object$rank <- object$qr$rank
    pred=predict.mlm(object, newdata = newdata, ...)
    if(ncol(pred)==1) pred=t(pred)
    pred
}

vcov.FRBmultireg <- function (object, ...) 
{
    if(!is.null(object$cov)) object$cov
    else cat("Not yet implemented\n")	 
}


