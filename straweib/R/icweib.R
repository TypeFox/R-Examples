icweib <-
function(L, R, data, strata="ALL", covariates=NULL) {
  Call <- match.call()
  if (is.name(Call$strata) | is.call(Call$strata)) {
  	straname <- deparse(Call$strata)
  } else {
  	straname <- "strata"
  }
  ## Get strata, left time, right time
  strata <- eval(substitute(strata), data, parent.frame())	
  Lt <- eval(substitute(L), data, parent.frame())	
  Rt <- eval(substitute(R), data, parent.frame())
  if (length(strata)==1) strata <- rep(strata, length(Lt))
  Lt <- as.vector(Lt)
  Rt <- as.vector(Rt)
  if (!is.numeric(Lt) | !is.numeric(Rt)) stop("L and R must be numeric variables")
  strata <- as.vector(strata)	
  
  ## Model Frame
  if (is.null(covariates)) {
  	mframe <- NULL
  	design <- NULL
  	beta.nm <- NULL
  	nbeta <- 0
  } else {
    mframe <-  model.frame(covariates, data=data, na.action=na.pass) 	
  }

  ## Compare lengths
  if (length(unique(c(length(Lt), length(Rt), length(strata), dim(mframe)[1])))!=1)
     stop("L, R, strata, and data must have same length")
     
  ## Remove missing
  nmiss <- complete.cases(cbind(Lt, Rt, strata, mframe))
  validtime <- (Lt >= 0) & (Rt > 0) & (Rt - Lt >= 0)
  validtime[is.na(validtime)] <- F
  keep <- nmiss & validtime
  delete <- which(!keep)
  ndel <- length(delete)	
  n <- sum(keep)
  Lt <- Lt[keep]
  Rt <- Rt[keep]
  strata <- strata[keep]
  mframe <- mframe[keep, , drop=F]
  
  ## Design matrix
  if (!is.null(covariates)) {
  	design <- model.matrix(covariates, data=mframe)[, -1, drop=F]
  	beta.nm <- colnames(design)
    nbeta <- dim(design)[2]
  }
	
  ## Design matrix for strata
  levels <- sort(unique(strata))
  nstr <- length(levels)
  dstrata <- outer(strata, levels, "==")*1
  stra.nm <- paste(c(rep("v", nstr), rep("u", nstr)), rep(levels, 2), sep=":")
  
  ## Time information
  allt <- c(Lt, Rt)
  maxT <- max(allt[allt!=Inf])
  meanT <- mean(allt[allt!=0 & allt!=Inf])
  vini <- -log(meanT)

  ## Get event status to see if any event
  et <- Lt==Rt
  event <- sum(et) > 0
  if (event) {
  	ie <- which(et)
  	dstratae <- dstrata[et, ]
  	designe <- design[et, ]
  	Rt[et] <- Inf  	
  }
  	
  ## Log-likelihood function
  LLt <- log(Lt)
  RRt <- log(Rt)	
  loglik <- function(parms) {
    v <- dstrata%*%parms[1:nstr]
    ev <- exp(v)
    u <- dstrata%*%parms[(nstr+1):(2*nstr)]
    if (nbeta > 0) {prog <- design%*%parms[-(1:(2*nstr))]} else {prog <- rep(0, n)}
    lik <- sum(log(exp(-exp(u + prog + ev*LLt)) - exp(-exp(u + prog + ev*RRt))))
    if (event) lik <- lik + sum(u[ie] + prog[ie] + v[ie] + (ev[ie]-1)*LLt[ie])
    return(-lik)	
  }
  	
  ## Gradient function		
  gradlik <-function(parms) {	
    v <- dstrata%*%parms[1:nstr]
    ev <- c(exp(v))
    u <- dstrata%*%parms[(nstr+1):(2*nstr)]
    if (nbeta > 0) {prog <- design%*%parms[-(1:(2*nstr))]} else {prog <- 0}
    SL <- exp(-exp(u + prog + ev*LLt))
    SR <- exp(-exp(u + prog + ev*RRt))
    SLL <- c(SL*exp(u + prog + ev*LLt))
    SRR <- c(ifelse(Rt==Inf, 0, SR*exp(u + prog + ev*RRt)))
    Dev <- colSums((cbind(ifelse(Lt==0, 0, LLt)*ev*dstrata, dstrata, design)*SLL - 
           cbind(ifelse(Rt==Inf, 0, RRt)*ev*dstrata, dstrata, design)*SRR)/c(SL-SR))
    if (event) Dev <- Dev - colSums(cbind((1+LLt[ie]*ev[ie])*dstratae, dstratae, designe))
    return(Dev)
  }

  ## Optimization
  parmi <- c(rep(0, nstr), rep(vini, nstr), rep(0, nbeta))
  q <- optim(parmi, loglik, gradlik, method="BFGS", hessian=T)
  if (q$convergence!=0)	warning("Full model not converged")
	
  ## Covariance matrix
  covm <- solve(q$hessian)	
  
  ## Fit NULL model, no covariates
  loglik0 <- function(parms) {
    v <- dstrata%*%parms[1:nstr]
    ev <- exp(v)
    u <- dstrata%*%parms[(nstr+1):(2*nstr)]
    lik <- sum(log(exp(-exp(u + ev*LLt)) - exp(-exp(u + ev*RRt))))
    if (event) lik <- lik + sum(u[ie] + v[ie] + (ev[ie]-1)*LLt[ie])
    return(-lik) 	
  }
  gradlik0 <- function(parms) {
    v <- dstrata%*%parms[1:nstr]
    ev <- c(exp(v))
    u <- dstrata%*%parms[(nstr+1):(2*nstr)]
    SL <- exp(-exp(u + ev*LLt))
    SR <- exp(-exp(u + ev*RRt))
    SLL <- c(SL*exp(u + ev*LLt))
    SRR <- c(ifelse(Rt==Inf, 0, SR*exp(u + ev*RRt)))
    Dev <- colSums((cbind(ifelse(Lt==0, 0, LLt)*ev*dstrata, dstrata)*SLL - 
           cbind(ifelse(Rt==Inf, 0, RRt)*ev*dstrata, dstrata)*SRR)/c(SL-SR))
    if (event) Dev <- Dev - colSums(cbind((1+LLt[ie]*ev[ie])*dstratae, dstratae))
    return(Dev)  
  }
  parmi0 <- c(rep(0, nstr), rep(vini, nstr))
  q0 <- optim(parmi0, loglik0, gradlik0, method="BFGS")
  	
  ## Fit the reduced model: shape parameters are all equal if there are more than 1 strata
  if (nstr > 1) {
    loglik1 <- function(parms) {
      v <- parms[1]
      ev <- exp(v)
      u <- dstrata%*%parms[2:(nstr+1)]
      if (nbeta > 0) {prog <- design%*%parms[-(1:(nstr+1))]} else {prog <- rep(0, n)}
      lik <- sum(log(exp(-exp(u + prog + ev*LLt)) - exp(-exp(u + prog + ev*RRt))))
      if (event) lik <- lik + sum(u[ie] + prog[ie] + v + (ev-1)*LLt[ie])
      return(-lik)	
	}	
	gradlik1 <-function(parms) {
      v <- parms[1]
      ev <- exp(v)
      u <- dstrata%*%parms[2:(nstr+1)]
      if (nbeta > 0) {prog <- design%*%parms[-(1:(nstr+1))]} else {prog <- 0}
      SL <- exp(-exp(u + prog + ev*LLt))
      SR <- exp(-exp(u + prog + ev*RRt))
      SLL <- c(SL*exp(u + prog + ev*LLt))
      SRR <- c(ifelse(Rt==Inf, 0, SR*exp(u + prog + ev*RRt)))
      Dev <- colSums((cbind(ifelse(Lt==0, 0, LLt)*ev, dstrata, design)*SLL - 
              cbind(ifelse(Rt==Inf, 0, RRt)*ev, dstrata, design)*SRR)/c(SL-SR))
      if (event) Dev <- Dev - colSums(cbind(1+LLt[ie]*ev, dstratae, designe))
      return(Dev)
	}
    parmi1 <- c(0, rep(vini, nstr), rep(0, nbeta))
    q1 <- optim(parmi1, loglik1, gradlik1, method="BFGS")
    if (q1$convergence!=0)	warning("Reduced model not converged")
	
	## Get likelihood ratio test statistic
    df <- nstr - 1
    testlik <- 2*(q1$value - q$value)
    plik <- 1 - pchisq(testlik, df)
    likratio <- data.frame(test="Likelihood Ratio", TestStat=testlik, df, p.value=plik)

	## Wald test statistic
    covv <- covm[1:nstr,1:nstr]
    estv <- q$par[1:nstr]
    testM <- t(sapply(1:(nstr-1), function(x) c(rep(0, x-1), c(1, -1), rep(0, nstr-x-1))))
    testwald <- t(testM%*%estv)%*%solve(testM%*%covv%*%t(testM))%*%(testM%*%estv)
    pwald <- 1 - pchisq(testwald, df)
    tWald <- data.frame(test="Wald", TestStat=testwald, df, p.value=pwald)
	
	## Combine test statistics
    stratatest <- rbind(tWald, likratio)		
    } else {
      stratatest <- NA
      q1 <- list(value=NA)
    }

	## report results
    rownames(covm) <- colnames(covm) <- c(stra.nm, beta.nm)
    logliks <- c(-q$value, -q1$value, -q0$value)
    names(logliks) <- c("full", "reduced", "null")
    if (nbeta > 0) {
	  beta.fit <- q$par[-(1:(2*nstr))]
	  beta.sd <- sqrt(diag(covm)[-(1:(2*nstr))])
	  beta.z <- beta.fit/beta.sd
	  p.value <- 2 * (1 - pnorm(abs(beta.z)))
	  coef <- data.frame(coefficient = beta.fit, SE = beta.sd, z = beta.z, p.value = p.value) 
	  rownames(coef) <- beta.nm    	
    } else {
      coef <- NA
    }
	weib <- data.frame(straname=straname, strata=levels, gamma = exp(q$par[1:nstr]),
	             lambda = exp(q$par[(nstr+1):(2*nstr)]), stringsAsFactors=F)
	ns <- c(n, nstr, nbeta, ndel)
	names(ns) <- c("nused", "nstrata", "ncovariates", "ndeleted")
	z <- list(loglik=logliks, coef=coef, weib=weib, stratatest=stratatest, cov=covm, ns=ns, delete=delete, maxT=maxT, q=q)
	class(z) <- "icweib"
	return(z)
}
