#' alr EM-based Imputation for Rounded Zeros
#' 
#' A modified EM alr-algorithm for replacing rounded zeros in compositional
#' data sets.
#' 
#' Statistical analysis of compositional data including zeros runs into
#' problems, because log-ratios cannot be applied.  Usually, rounded zeros are
#' considerer as missing not at random missing values. The algorithm first
#' applies an additive log-ratio transformation to the compositions. Then the
#' rounded zeros are imputed using a modified EM algorithm.
#' 
#' @param x Compositional data
#' @param pos Position of the rationing variable for alr transformation
#' @param dl Detection limit for each part
#' @param eps convergence criteria
#' @param maxit maximum number of iterations
#' @param bruteforce if TRUE, imputations over dl are set to dl. If FALSE,
#' truncated (Tobit) regression is applied.
#' @param method either \dQuote{lm} (default) or \dQuote{MM}
#' @param step if TRUE, a stepwise (AIC) procedure is applied when fitting
#' models
#' @param nComp if determined, it fixes the number of pls components. If
#' \dQuote{boot}, the number of pls components are estimated using a
#' bootstraped cross validation approach.
#' @param R number of bootstrap samples for the determination of pls
#' components. Only important for method \dQuote{pls}.
#' @param verbose additional print output during calculations.
#' @return \item{xOrig }{Original data frame or matrix} \item{xImp }{Imputed
#' data} \item{wind }{Index of the missing values in the data} \item{iter
#' }{Number of iterations} \item{eps }{eps}
#' @author Matthias Templ and Karel Hron
#' @seealso \code{\link{impRZilr}}
#' @keywords manip multivariate
#' @export
#' @importFrom MASS stepAIC
#' @importFrom robustbase lmrob
#' @examples
#' 
#' data(arcticLake)
#' x <- arcticLake
#' ## generate rounded zeros artificially:
#' x[x[,1] < 5, 1] <- 0
#' x[x[,2] < 47, 2] <- 0
#' xia <- impRZalr(x, pos=3, dl=c(5,47), eps=0.05)
#' xia$xImp
#' 
impRZalr <- function(x, pos=ncol(x), dl=rep(0.05, ncol(x)-1), 
                     eps=0.0001, maxit=50, bruteforce=FALSE, 
                     method="lm", step=FALSE, nComp = "boot", R=10,
                     verbose=FALSE) {

#   x <- xorig <- genVarsSimple(n=50,p=30)
#   qu <- apply(x, 2, quantile, 0.05)
#   for(i in 1:(ncol(x)-2)){
#      x[x[,i] < qu[i], i] <- 0
#   }
#   pos <- ncol(x)-1
#   x <- data.frame(x)
#   dl <- qu
#   method="lm"
#   step=FALSE
#   bruteforce=FALSE
  
  ## some checks:	
  if(is.matrix(x)) stop("convert to data.frame first")
  stopifnot(all(x[,pos] != 0 & length(which(is.na(x[,pos]))) == 0))
  if(!any(is.na(x)) && !any(x==0) ) stop("No zeros or missing values in the data")
  if(method=="pls" & ncol(x)<5) stop("too less variables/parts for method pls")
  if(is.null(nComp)){
    pre <- FALSE
    nC <- NULL
  } else if(nComp=="boot"){
    nC <- integer(ncol(x))
    pre <- TRUE
  } else if(length(nComp) == ncol(x)){
    nC <- nComp
    pre <- FALSE
  } else  {
    pre <- FALSE	
  }
  
  ## zeros to NA:
  x[x==0] <- NA
  
  ## transformation:
  xa <- addLR(x, ivar=pos)
  xax <- xa$x
  w <- is.na(xa$x)
  
  ## dl --> phi
  m <- matrix(rep(dl, each=nrow(x)), ncol=length(dl))
  phi <- log(m/x[,pos]) 
  phi <- phi[,-pos]
  xOrig <- x

  it <- 0
  d <- 99999999
  
  ## initialisation:
  xax[w] <- phi[w] * 2 / 3
  
  ## start the EM:
  it <- 0
  amountMiss <- length(which(w))
  while( d > eps & it <= maxit ){
	it <- it + 1
	yold <- xax
  n2 <- nrow(xax)-ncol(xax)+1
  cn <- colnames(x)
  n <- nrow(x)
  for(i in 1:ncol(xax)){ 
    response <- xax[,i]
    predictors <- as.matrix(xax[,-i])
		if(method=="lm" && !step){
	    lm1 <- lm(response ~ predictors)
  			yhat <- predict(lm1, new.data=predictors)	  
  			s <- sd(lm1$res, na.rm=TRUE)
	#  		s <- sqrt(sum(lm1$res^2)) / n2	
	  }
    if(method=="pls" && !step){    
      if(it == 1 & pre){ ## evaluate ncomp.
        nC[i] <- bootnComp(predictors,
                           y=response, R, 
                           plotting=FALSE)$res #$res2
      }
      if(verbose) cat("\n ncomp for part",i,":",nC[i])
      reg1 <- mvr(as.matrix(response) ~ as.matrix(predictors), ncomp=nC[i],
                  method="simpls")
      yhat <- predict(reg1, new.data=data.frame(predictors), ncomp=nC[i])
      s <- sqrt(sum(reg1$res^2)/n) 
      
      lm1 <- lm(response ~ predictors)
      yhat <- predict(lm1, new.data=predictors)	  
    		s <- sd(lm1$res, na.rm=TRUE)
      # s <- sqrt(sum(lm1$res^2)) / n2	
    }
    if(method=="lm" && step){
	  	lm1 <- lm(response ~ predictors) 
	    lm1 <- MASS::stepAIC(lm1, trace=FALSE)
	  	yhat <- predict(lm1, new.data=predictors)	  
    	s <- sd(lm1$res, na.rm=TRUE)
  	#	s <- sqrt(sum(lm1$res^2)) / n2			
	 	}
    if(method=="MM" && !step) {
  		lm1 <- robustbase::lmrob(response ~ predictors)
  		yhat <- predict(lm1, new.data=predictors)	
  		if(any(is.na(yhat))) stop("NA in yhat")
  		if(any(yhat=="Inf")) stop("Inf in yhat")	
  		s <- lm1$s		
  	}
	  if(method=="MM" && step) {
	    stop("robust + stepwise is not implemented until now.")
	   }	
	  ex <- (phi[,i] - yhat)/s 
		tobit <- s*dnorm(ex)/pnorm(ex)
		tobit <- ifelse(tobit=="NaN", 0, tobit)
		tobit <- ifelse(is.na(tobit), 0, tobit)	
		tobit <- ifelse(tobit =="Inf", 0, tobit)	
		tobit <- ifelse(tobit =="-Inf", 0, tobit)	
		yhat2 <- yhat - tobit
	  # check if we are under the DL:
    if(any(yhat2[w[,i]] >= phi[w[,i],i])) stop(paste("values above the DL are imputed. \n column",i))
	   if(bruteforce){
          xax[w[,i],i] <- ifelse(yhat2[w[,i]] >= phi[w[,i],i], phi[w[,i],i], yhat2[w[,i]]) 
	   } else {
       xax[w[,i],i] <- yhat2[w[,i]] 
	   }
	}
	d <- sum(abs(xax - yold))/amountMiss
}

  ## backtransform:	
  xa$x.alr <- xax
  ximp <- suppressWarnings(addLRinv(xa)) 
  
  ## result:
  res <- list(x=ximp, wind=NULL, iter=it, eps=eps) 
  invisible(res)
}


