phm <- function(data, ...)UseMethod("phm")

phm.default <- function(data = NULL, subset=NULL, 
                        TP="TP", FN="FN", FP="FP", TN="TN",  
                        correction = 0.5, correction.control = "all", 
                        hetero = TRUE, estimator = "APMLE", l = 100, ...){  
  if(estimator == "APMLE"){estimator <- "Adjusted Profile Maximum Likelihood"}
  DNAME <- deparse(substitute(x))
 
  
  stopifnot(is.numeric(correction), 0 <= correction,  
            correction.control %in% c("all", "single", "none"),
            is.numeric(TP) | (is.character(TP) & length(TP) == 1),
            is.numeric(FP) | (is.character(FP) & length(FP) == 1),
            is.numeric(TN) | (is.character(TN) & length(TN) == 1),
            is.numeric(FN) | (is.character(FN) & length(FN) == 1))
  
  if(!is.null(data) & is.character(c(TP,FP,TN,FN))){
    X <- as.data.frame(data)
    origdata <- data
    TP <- getElement(X,TP)
    FN <- getElement(X,FN)
    FP <- getElement(X,FP)
    TN <- getElement(X,TN)
  }
  
  if(is.null(data) & !is.character(c(TP,FP,TN,FN))){
    origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
  }
  
  freqdata <- cbind(TP,FN,FP,TN)
  checkdata(freqdata)
  
  k <- length(TP)  
  
  ## apply continuity correction to _all_ studies if one contains zero
  if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0))
  {TP <- TP + correction;
   FN <- FN + correction;
   FP <- FP + correction;
   TN <- TN + correction}}
  if(correction.control == "single"){
    correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0)) * 
      correction
    TP <- correction + TP
    FN <- correction + FN
    FP <- correction + FP
    TN <- correction + TN
  }

  theta_A=numeric(l)
	tau_sq_A=numeric(l)

  if(hetero){

	y <- TP
	m <- TP+FN
	x <- FP
	n <- FP+TN
	z <- log(y)-log(m)
	s2 <- 1/y-1/m
	w <- log(x)-log(n)
	t2 <- 1/x-1/n

  #initialize
	theta_A[1] <- 0.5
	tau_sq_A[1]<-0
	theta_hat <- z/w
  
	for (j in 2:l){
		sigma_sq <- t2*theta_A[j-1]^2+s2
		v <- 1/(sigma_sq/w^2+tau_sq_A[j-1])
		theta_Aj <-rep(theta_A[j-1],k)
		theta_A[j]<- sum(theta_hat*v)/sum(v-(theta_hat-theta_Aj)^2*v^2*(t2/w^2)+(t2/w^2)*v)
		tau_sq_A[j] <- sum(((theta_hat-theta_A[j-1])^2-sigma_sq/w^2)*v^2)/sum(v^2)	

    if(l == j){warning("Reached maximum number of iterations!")}  
    if(abs(theta_A[j]-theta_A[j-1])<0.00000001)break
	}
		sigma_sq <- t2*theta_A[j]^2+s2
			if(tau_sq_A[j]<0){tau_sq_A[j]<-0}
		v <- 1/(sigma_sq/w^2+tau_sq_A[j])
		theta_Aj <-rep(theta_A[j],k)
		a <- (theta_hat-theta_Aj)*v+(theta_hat-theta_Aj)^2*v^2*(t2/w^2)*theta_Aj-(t2/w^2)*v*theta_Aj
		c  <- 0.5*(theta_hat-theta_Aj)^2*v^2-0.5*v
		D <- sum(a^2)*sum(c^2)-(sum(a*c))^2
		var_theta_A <- sum(c^2)/D

		var_tau_sq_A <- sum(a^2)/D

    max_log_theta_A <-  -0.5*sum(log(sigma_sq/(w^2) + tau_sq_A[j] )) - 
      (sum(0.5*((theta_hat - theta_A[j])^2/(sigma_sq/(w^2) + tau_sq_A[j]))))
		score_A <- sum(a)
		score_AT <- sum(c)
		chi_sq_A <- sum((theta_hat-theta_Aj)^2*v)
} # End hetero
  
if(!hetero){
    y <- TP
	  m <- TP+FN
	  x <- FP
	  n <- FP+TN
		z <- log(y)-log(m)
		s2 <- 1/y-1/m
		w <- log(x)-log(n)
		t2 <- 1/x-1/n

		theta_A[1] <- sum(z*w)/sum(w^2)
		for (j in 2:l)
		{
		sigma_sq <- t2*theta_A[j-1]^2+s2
		v <- 1/sigma_sq
		theta_A[j]<- sum(z*w*v)/sum(w^2*v-(z-w*theta_A[j-1])^2*v^2*t2+t2*v) 
	if(abs(theta_A[j]-theta_A[j-1])<0.00000001)
		break
		}
		sigma_sq <- t2*theta_A[j]^2+s2
		v <- 1/sigma_sq

		a <- (z-w*theta_A[j])*w*v+(z-w*theta_A[j])^2*v^2*t2*theta_A[j]-t2*v*theta_A[j]

		var_theta_A <- 1/sum(a^2)

		max_log_theta_A <-  -0.5*sum(log(sigma_sq)) -(sum(0.5*((z-w*theta_A[j])^2/sigma_sq)))

		score_A <- sum(a)
		chi_sq_A <- sum((z-w*theta_A[j])^2*v)
    
    var_tau_sq_A <- NULL
    score_AT <-NULL
} # End !hetero
  
  
  if(hetero){
  coefficients <- c(theta_A[j], tau_sq_A[j])
  names(coefficients) <- c("theta", "taus_sq")
  }else{
  coefficients <- c(theta_A[j])
  names(coefficients) <- c("theta")
  }
  
  if(hetero){
  vcov <- diag(c(var_theta_A, var_tau_sq_A))
  colnames(vcov) <- c("theta", "taus_sq")
  rownames(vcov) <- c("theta", "taus_sq")
  }else{
  vcov <- matrix(var_theta_A, nrow = 1, ncol = 1)
  colnames(vcov) <- "theta"
  rownames(vcov) <- "theta"
  }
  
  score <- c(score_A, score_AT)
  if(hetero){
  names(score) <-c("score_theta", "score_tau_sq")
  }else{
  names(score) <-c("score_theta")    
  }
  
  parameter <- (1+hetero)
  names(parameter) <- "df"
  names(chi_sq_A) <- "Chi-square"
  chi_sq_test <- list(statistic = chi_sq_A, parameter = parameter,
                      method =paste(c("Chi-square goodness of fit test (", 
                                estimator,ifelse(hetero, " under heterogeneity)", 
                                                 " under homogeneity)")), collapse = "", sep=""),
                      data.name = DNAME, 
                      p.value = 1-pchisq(chi_sq_A,k-(1+hetero)))
  class(chi_sq_test) <- "htest"
  
  ll <- max_log_theta_A
  attr(ll, "df") <- (1+hetero)
  
  output <- list(estimator = estimator, hetero = hetero, coefficients = coefficients,
                 vcov = vcov, logLik = ll, chi_sq_test = chi_sq_test,
                 iterations = j, call = match.call(), nobs = k, 
                 data = origdata)
  class(output) <- "phm"
  return(output)  
}

summary.phm <- function(object, level = 0.95, ...)
{
output <-list(object = object, level = level, AUC = AUC(object, level = level))
class(output) <- "summary.phm"
return(output)
}

print.summary.phm <- function(x, ...)
{
cat("Call:\n")
print(x$object$call)
cat("\n")
print(cbind(Estimate = coef(x$object), confint(x$object, level = x$level)))
cat("\n")
cat(c("Log-likelihood:", round(logLik(x$object),3), 
      "on",attr(logLik(x$object), "df"), "degrees of freedom\n"))
cat(c("AIC: ", round(AIC(x$object),1), "\n"))
cat(c("BIC: ", round(AIC(x$object, k = log(x$object$nobs)),1), "\n"))

print(x$object$chi_sq_test)

cat("\n")
auc <- c(x$AUC$AUC, x$AUC$ci, x$AUC$pAUC, x$AUC$pci)
names(auc)[1] <- "AUC"
names(auc)[4] <- "pAUC"
print(round(auc,3))
}

logLik.phm <- function(object, ...){object$logLik}
vcov.phm <- function(object, ...){object$vcov}


print.phm <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

sroc.phm <- function(fit, fpr = 1:99/100, ...){
sens <-  fpr^coef(fit)[1]
return(cbind(fpr = fpr, sens = sens))
}

srocphm <- function(theta, fpr = 1:99/100, ...){
sens <-  fpr^theta
return(cbind(fpr = fpr, sens = sens))
}


plot.phm <- function(x, extrapolate = FALSE, confband = TRUE, level = 0.95, 
                     ylim = c(0,1), xlim = c(0,1),
                     sroclty = 1, sroclwd = 1,
                     confbandlty = 2, confbandlwd = 0.5, ...){
  FP <- x$data$FP
  negatives <- FP + x$data$TN
  FPR <- FP/negatives
  
  if(extrapolate){bound = c(0,1)}
  if(!extrapolate){bound = c(min(FPR), max(FPR))}
  plot(c(2,2), ylim = ylim, xlim = xlim, 
       xlab = "False Positive Rate", ylab = "Sensitivity")

  srocmat <- sroc(x)
  lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",], 
        lty = sroclty, lwd = sroclwd, ...)
  if(confband){conftheta <- confint(x)[1,]
               if(conftheta[1]>0){
                 srocmat <- srocphm(conftheta[1])
                 lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",], 
                      lty = confbandlty, lwd = confbandlwd, ...)}
                 srocmat <- srocphm(conftheta[2])
                 lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",], 
                      lty = confbandlty, lwd = confbandlwd, ...)
               }
}