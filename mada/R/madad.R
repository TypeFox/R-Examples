# function, output of class "madad"
madad <- function(x=NULL, TP, FN, FP, TN, level = 0.95, correction = 0.5, 
                         correction.control = "all", method = "wilson", 
                  yates = TRUE, suppress = TRUE, ...){
  if(is.data.frame(x)){names <- x$names}else{names <- NULL}
  alpha<-1-level
	kappa<-qnorm((1-alpha/2))
  DNAME <- deparse(substitute(x))
  if(!is.null(x)){
    X <- as.data.frame(x)
    TP <- X$TP
    FN <- X$FN
    FP <- X$FP
    TN <- X$TN
    origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
  }
  if(is.null(x)){origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)}

  checkdata(origdata)

  ## apply continuity correction to _all_ studies if one contains zero
  if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0)){TP <- TP+correction;
    					FN <- FN + correction;
							FP <- FP + correction;
							TN <- TN + correction}}
  if(correction.control == "single"){
	  correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0))*correction
    TP <- correction + TP
	  FN <- correction + FN
	  FP <- correction + FP
	  TN <- correction + TN
	}
  
  
  number.of.pos<-TP+FN
  number.of.neg<-FP+TN
  sens.ci<-binomCIvector(TP,number.of.pos, conf.level  = level, method = method)
  colnames(sens.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))
  spec.ci<-binomCIvector(TN,number.of.neg, conf.level  = level, method = method)
  colnames(spec.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))
  fpr.ci<-binomCIvector(FP,number.of.neg, conf.level  = level, method = method)
  colnames(fpr.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))


  sens <- TP/number.of.pos
  spec <- TN/number.of.neg
  fpr <- FP/number.of.neg

  cor_da <- cor(as.matrix(data.frame(sens = sens, spec = spec, fpr = fpr)))
  
	posLR<-sens/fpr
	negLR<-(1-sens)/(1-fpr)
	
	se.lnposLR <- sqrt(1/TP + 1/FP - 1/number.of.pos - 1/number.of.neg)
	se.lnnegLR <- sqrt(1/TN + 1/FN - 1/number.of.pos - 1/number.of.neg)
	posLR.ci<-cbind(exp(-kappa*se.lnposLR),exp(kappa*se.lnposLR))
	posLR.ci<-posLR*posLR.ci
  colnames(posLR.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))

	negLR.ci<-cbind(exp(-kappa*se.lnnegLR),exp(kappa*se.lnnegLR))
	negLR.ci<-negLR*negLR.ci
	colnames(negLR.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))

  DOR<-posLR/negLR
  se.lnDOR<-sqrt(1/TP + 1/TN + 1/FN + 1/FP)
	DOR.ci<-cbind(exp(-kappa*se.lnDOR),exp(kappa*se.lnDOR))
	DOR.ci<-DOR*DOR.ci
  colnames(DOR.ci) <- c(paste(100*alpha/2, "%", sep="", collapse = ""), paste(100*(1-alpha/2), "%", sep="", collapse = ""))
  
  
  sens.htest <- switch(suppress, 
                       suppressWarnings(prop.test(TP, number.of.pos, 
                                                  correct = yates)),
                       prop.test(TP, number.of.pos, correct = yates))
                       
                       
  spec.htest <- switch(suppress, 
                       suppressWarnings(prop.test(TN, number.of.neg, 
                                                  correct = yates)),
                       prop.test(TN, number.of.neg, 
                                 correct = yates))
  
  output <- list(sens = list(sens = sens, sens.ci = sens.ci),
                 spec = list(spec = spec, spec.ci = spec.ci),
                 fpr = list(fpr = fpr, fpr.ci = fpr.ci),
                 sens.htest = sens.htest,
                 spec.htest = spec.htest,
                 posLR = list(posLR = posLR, posLR.ci = posLR.ci, se.lnposLR = se.lnposLR),
                 negLR = list(negLR = negLR, negLR.ci = negLR.ci, se.lnnegLR = se.lnnegLR),
                 DOR = list(DOR = DOR, DOR.ci = DOR.ci, se.lnDOR = se.lnDOR),
                 cor_sens_fpr = cor(sens,fpr),
                 level = level, method = method, names = names,
                 nobs = nrow(origdata), data = origdata, data.name = DNAME,
                 correction = correction, correction.control = correction.control)
  class(output) <- "madad"
  output
}

print.madad <- function(x, digits = 3, ...){
  cat("Descriptive summary of", x$data.name, "with", x$nobs, "primary studies.\n")
  cat("Confidence level for all calculations set to", 100*x$level, "%\n")
  cat("Using a continuity correction of", x$correction, "if applicable \n")
  cat("\n")
  
  
  cat("Diagnostic accuracies \n")
  output1a <- round(cbind(x$sens$sens, x$sens$sens.ci, x$spec$spec, x$spec$spec.ci), digits)
  rownames(output1a) <- x$names
  colnames(output1a)[c(1,4)] <- c("sens", "spec")
  print(output1a)
  
  format.test.result<-function(x)
  {
	out <- character()
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 
            4))), ", ")
    if (!is.null(x$parameter)) 
        out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 
            3))), ", ")
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        out <- c(out, paste("p-value =", fp))}
    return(out)
	} # End of function format.test.result

  sens.result<-format.test.result(x$sens.htest)
	spec.result<-format.test.result(x$spec.htest)
  
	output1b<-paste(c("\nTest for equality of sensitivities: \n",sens.result, "\n",
                   "Test for equality of specificities: \n",spec.result, "\n\n"), collapse="")
	cat(output1b)
  
  cat("\n")
  cat("Diagnostic OR and likelihood ratios \n")
  output2 <- round(cbind(x$DOR$DOR, x$DOR$DOR.ci, x$posLR$posLR, x$posLR$posLR.ci,
                         x$negLR$negLR, x$negLR$negLR.ci), digits)
  rownames(output2) <- x$names
  colnames(output2)[c(1,4,7)] <- c("DOR", "posLR", "negLR")
  print(output2)
  
  cat("\n")
  cat("Correlation of sensitivities and false positive rates: \n")
  print(round(CIrho(x$cor_sens_fpr, x$nobs), digits)[1,])
  return(invisible(NULL))
  }
