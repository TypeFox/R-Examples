summary.glmssn <-
function(object, ...)
{
  catCall =  object$args$call 

	effnames <- object$sampinfo$effnames
	setzero <- object$sampinfo$setzero
	setNA <- object$sampinfo$setNA
	if(length(setNA)==1) setNA <- FALSE
	setNA2 <- object$sampinfo$setNA2
	if(length(setNA2)==1) setNA2 <- FALSE
	b.hat <- object$estimates$betahat
	bhat.se <- sqrt(diag(object$estimates$covb))
	n.allxy <- object$sampinfo$obs.sample.size
	p <- object$sampinfo$rankX

  if(any(rownames(b.hat) %in% effnames == FALSE)) {
            ## dataXY issue
  stop(cat("glmssn has computed estimates for",rownames(b.hat),"but the summary command expects estimates for",effnames,collapse=" "))
        }

	bhat.0.NA <- rep(NA, times = length(effnames))
	bhat.0.NA[setzero] <- 0
	bhat.0.NA[!setzero & !setNA2] <- b.hat
	NAvec <- rep(NA, times = length(effnames))
	bhatse.0.NA <- NAvec
	bhatse.0.NA[!setzero & !setNA2] <- bhat.se
	tvec <- NAvec
	tvec[!setzero & !setNA2] <- b.hat/bhat.se
	pvec <- NAvec
	pvec[!setzero & !setNA2] <-
			round(100000*(1 - pt(abs(b.hat/bhat.se), df = n.allxy - p))*2)/100000
	fixed.eff.est <- data.frame(FactorLevel = effnames, Estimate = bhat.0.NA,
			std.err = bhatse.0.NA, t.value = tvec, prob.t = pvec)
	fixed.effects.estimates = fixed.eff.est

  if(object$args$algorithm=="orig") {
    covmodels <- object$estimates$theta
    covmodels <- data.frame(Covariance.Model=attributes(
		  covmodels)$terms,Parameter=attributes(covmodels)$type,
			Estimate=covmodels)
  }

  if(object$args$algorithm=="multi") {
    covmodels <- object$estimates$theta
    covmodels <- data.frame(Covariance.Model=attributes(covmodels)$terms,
      Parameter=attributes(covmodels)$type,Estimate=covmodels)
  }

  if(object$args$algorithm=="regress") {
    cat("improve this summary\n")
    covmodels <- object$estimates$theta
  }
	 
  res = getSSNdata.frame(residuals(object))[,"_resid_"]

  Warnlog <- object$estimates$Warnlog

  outpt = list(catCall = catCall,
    fixed.effects.estimates = fixed.effects.estimates,
    covariance.parameter.estimates = covmodels,
    res = res,
    Rsquared = GR2(object),
		Warnings = Warnlog)
  class(outpt) <- "summary.glmssn"
  outpt
}

print.glmssn <- function(x,...) {
    print(summary(x,...))
}

print.summary.glmssn <- function(x, 
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$catCall), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
  cat("\nResiduals:\n")
  resQ = c(min(x$res), quantile(x$res, p = c(.25, .5, .75), 
    na.rm = TRUE), max(x$res))
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)
  cat("\nCoefficients:\n")
  coef1 = x$fixed.effects.estimates
  coefs = coef1[,2:5]
  row.names(coefs) = coef1[,1]
  colnames(coefs) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
  cat("\nCovariance Parameters:\n")
  cpe = x$covariance.parameter.estimates
  dig = -min(0, -2 + as.numeric(strsplit(formatC(min(cpe[,3]), digits = 3, 
    format = "E"),"E")[[1]][2]))
  print(data.frame(cpe[1:2], Estimate = formatC(cpe[,3],
    digits = dig, format = "f", flag = "0")), row.names = FALSE)
  rse = sqrt(sum(cpe[cpe[,"Parameter"] == "parsill","Estimate"]))
  cat("\nResidual standard error:", rse)
  cat("\nGeneralized R-squared:", x$Rsquared)

  cat("\n")
  invisible(x)
}
