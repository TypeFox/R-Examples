
##########################################################################
###########  2 STATES
##########################################################################


######################################################################
###########       DISTRIBUTION CHOICE
######################################################################
########          RETURNS THE DENSITY OF THE CHOSEN DISTRIBUTION
##### arguments: dist=character indicating any 2-parameter distribution
#####            implemented in R
##### value: the chosen density function
######################################################################
dens <- function(dist){
  fun <- get(paste0("d", dist))
  return(fun)
}


######################################################################
########          RETURNS THE QUANTILES OF THE CHOSEN DISTRIBUTION
##### arguments: dist=character indicating any 2-parameter distribution
#####            implemented in R
##### value: the chosen quantile function
######################################################################
quant <- function(dist){
  fun <- get(paste0("q", dist))
  return(fun)
}


######################################################################
########          RETURNS THE DISTRIBUTION FUNCTION OF THE CHOSEN DISTRIBUTION
##### arguments: dist=character indicating any 2-parameter distribution
#####            implemented in R
##### value: the chosen distribution function
######################################################################
p <- function(dist){
  fun <- get(paste0("p", dist))
  return(fun)
}

######################################################################
########          GENERATES RANDOM NUMBERS OF THE CHOSEN DISTRIBUTION
##### arguments: dist=character indicating any 2-parameter distribution
#####            implemented in R
##### value: the chosen r function
######################################################################
rand <- function(dist){
  fun <- get(paste0("r", dist))
  return(fun)
}


############################################################################
########          BETA OR R FUNCTION (COST-MINIMISING SLOPE)
##### arguments:  rho=disease prevalence
#####		  	      costs=cost matrix
##### value: cost-minimising slope
############################################################################
slope <- function(rho, costs=matrix(c(0,0,1,(1-rho)/rho), 2, 2, byrow=TRUE)){
  c.t.pos <- costs[1, 1]
  c.t.neg <- costs[1, 2]
  c.f.pos <- costs[2, 1]
  c.f.neg <- costs[2, 2]
  beta <- ((1 - rho)/rho) * ((c.t.neg - c.f.pos)/(c.t.pos - c.f.neg))
  return(beta)
}


############################################################################
########          CHECKING POSITIVE SQUARE ROOT
##### arguments: k1=vector containing the healthy sample values 
#####		         k2=vector containing the diseased sample values
##### 	 	       rho=disease prevalence
#####		  	     costs=cost matrix
##### value: discriminant of the expression for the optimum threshold
############################################################################
sqroot <- function(k1, k2, rho, costs){
  ctrl <- (mean(k2) - mean(k1))^2 + 2*log((sd(k2)/sd(k1))*slope(rho,costs))*(var(k2) - var(k1))
  return(ctrl)
}


############################################################################
##### THRESHOLD DEPENDING ON THE UNDERLYING DISTRIBUTION OF THE POPULATIONS
##### COST-MINIMISING THRESHOLD ONE-VARIABLE EQUATION
##### DENSITY RATIO FORMULA
##### arguments:  p=vector containing the interval extremes at which the ?uniroot? function will look for 
#####                         the one-variable equation solution  
#####		          dist1=choose the healthy distribution
#####		   	      dist2=choose the diseased distribution			
#####             par1.1=healthy population mean
#####		   	      par1.2=healthy population standard deviation
#####             par2.1=diseased population mean
#####		          par2.2=diseased population standard deviation
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: the threshold one-variable equation evaluated at p
############################################################################
DensRatio2 <- function(p, dist1, dist2, par1.1, par1.2, par2.1, par2.2, rho, costs) {
  ratio <- (dens(dist2)(p,par2.1,par2.2)/dens(dist1)(p,par1.1,par1.2)) - slope(rho,costs)
  return(ratio)
}


############################################################################
##### "uniroot" function looks for the one-variable equation solution  
##### arguments: 	q1=probability of the "left" distribution in order to determine a "low" quantile.
##### 	 	        q2=probability of the "right" distribution in order to determine a "high" quantile.
#####			        dist1, dist2=choose the healthy and the diseased distributions, respectively,
#####             between the following 2-parameter distributions: "beta", "cauchy",
#####             "chisq" (chi-squared), "gamma", "lnorm" (lognormal), "logis" (logistic), "norm" (normal)
#####             and "weibull".
#####			        par1.1=first parameter of the distribution chosen for the healthy population
#####		   	      par1.2=second parameter of the distribution chosen for the healthy population
#####             par2.1=first parameter of the distribution chosen for the diseased population
#####		          par2.2=second parameter of the distribution chosen for the diseased population
#####		          rho=disease prevalence
#####		  	      costs=cost matrix
#####             tol=tolerance to be used in uniroot. Default, 10^(-8).
##### value: an object of class "thresTH2" containing the threshold one-variable equation solution,
#####        the prevalence, the cost matrix and the slope.      
############################################################################
thresTH2 <- function(dist1, dist2, par1.1, par1.2, par2.1, par2.2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), q1=0.05, q2=0.95, tol=10^(-8)){
  # error handling
  if (!(rho > 0 & rho < 1)){
    stop("The disease prevalence rho must be a number in (0,1)")
  }
  if (!is.matrix(costs)){
    stop("'costs' must be a matrix")
  }
  if (dim(costs)[1] != 2 | dim(costs)[2] != 2){
    stop("'costs' must be a 2x2 matrix")
  }
  # function
  costs.origin <- costs
  rho.origin <- rho
  median1 <- quant(dist1)(0.5, par1.1, par1.2)
  median2 <- quant(dist2)(0.5, par2.1, par2.2)
  if(median1 > median2){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- par2.1; par2.1 <- par1.1; par1.1 <- g
    f <- par2.2; par2.2 <- par1.2; par1.2 <- f
    auxdist <- dist2; dist2 <- dist1; dist1 <- auxdist
  }
  p1 <- quant(dist1)(q1, par1.1, par1.2)
  p2 <- quant(dist2)(q2, par2.1, par2.2)
  # threshold estimate
  cut.t <- uniroot(DensRatio2,c(p1,p2),tol=tol,dist1,dist2,par1.1,par1.2,par2.1,par2.2,rho,costs)$root
  # add slope
  beta <- slope(rho, costs)
  # results
  re <- list(thres=cut.t, prev=rho.origin, costs=costs.origin, R=beta, method="theoretical")
  # return
  class(re) <- "thresTH2"
  return(re)
}


############################################################################
# Print function for class "thresTH2"
############################################################################

print.thresTH2 <- function(x, ...){
    cat("\nThreshold:", x$thres)
    cat("\n")
    cat("\nParameters used")
    cat("\n  Disease prevalence:", x$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$costs)
    cat("\n  R:", x$R)
    cat("\n")
}


##########################################################################
##### Second derivative of the cost function (2 states): just a check
##### arguments:  x = thres2 object
##### value:      the second derivative function
##########################################################################
secondDer2 <- function(x){
  # error handling
  if (class(x) != "thres2"){
    stop("'x' must be of class 'thres2'")
  }
  # function
  if (x$T$method=="empirical"){
    stop("'x' has been computed with method='empirical', cannot compute the second derivative of the cost function")
  }else if (x$T$method=="equal" | x$T$method=="unequal"){
    k1 <- x$T$k1
    k2 <- x$T$k2
    rho <- x$T$prev
    costs <- x$T$costs
    Thr <- x$T$thres
    par1.1 <- mean(k1)
    par2.1 <- mean(k2)
    par1.2 <- sd(k1)
    par2.2 <- sd(k2)
    n1 <- length(k1)
    n2 <- length(k2)
    beta <- slope(rho,costs)
    de <- (n1+n2)*rho*(costs[1,1]-costs[2,2])
    der <- de/sqrt(2*pi)*(((Thr-par2.1)/par2.2^3)*exp(-(Thr-par2.1)^2/(2*par2.2^2))-beta*((Thr-par1.1)/par1.2^3)*exp(-(Thr-par1.1)^2/(2*par1.2^2)))
  }else if (x$T$method=="parametric"){
    dist1 <- x$T$dist1
    dist2 <- x$T$dist2
    pars1 <- x$T$pars1
    pars2 <- x$T$pars2
    densratio <- function(y){
      DR <- dens(dist2)(y, pars2[1], pars2[2])/dens(dist1)(y, pars1[1], pars1[2])-x$T$R
      return(DR)
    }
    der <- grad(densratio, x$T$thres)
  }
  return(der)
}


###########################################################################
############      ASSUMING EQUAL VARIANCES
##########################################################################
############      FUNCTION OF POOLED VARIANCE
##########################################################################
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
##### value:      the pooled variance of the two samples
##########################################################################
varPooled <- function(k1, k2){
  n1 <- length(k1)
  n2 <- length(k2)
  pool <- ((n1 - 1)*var(k1) + (n2 - 1)*var(k2))/(n1 + n2 - 2)
  return(pool)
}

##########################################################################
############   1) FUNCTION OF THRESHOLD - JUND (EQUAL VARIANCES)
##########################################################################
##### arguments:  k1=vector containing the healthy sample values 
#####		         	k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: the threshold estimated by the two samples considering the population variances equal
##########################################################################
thresEq2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE)){
  costs.origin <- costs
  k1.origin <- k1
  k2.origin <- k2
  rho.origin <- rho
  if(mean(k1)>mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  beta <- slope(rho, costs)
  # threshold estimate
  cut <- (2*varPooled(k1,k2)*log(beta) - (mean(k1)^2 - mean(k2)^2))/(2*(mean(k2) - mean(k1)))
  # returning results
  re <- list(thres=cut, prev=rho.origin, costs=costs.origin, R=beta, method="equal", k1=k1.origin, k2=k2.origin)
  return(re)
}



##########################################################################
###########    2) FUNCTION OF THRESHOLD - JUND (UNEQUAL VARIANCES)
##########################################################################
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		        	rho=disease prevalence
#####		        	costs=cost matrix
##### value: the threshold estimated by the two samples considering the population variances unequal
##########################################################################
thresUn2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE)){
  costs.origin <- costs
  k1.origin <- k1
  k2.origin <- k2
  rho.origin <- rho
  if(mean(k1)>mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  ctrl <- sqroot(k1, k2, rho, costs) 	
  if (ctrl>=0){
    # threshold estimate
    cut <- (var(k2)*mean(k1) - var(k1)*mean(k2) + sd(k1)*sd(k2)*sqrt(ctrl))/(var(k2) - var(k1))}
  else{
    cut <- NA
    warning("Negative discriminant; cannot solve the second-degree equation")
  }
  beta <- slope(rho, costs)
  # returning results
  re <- list(thres=cut, prev=rho.origin, costs=costs.origin, R=beta, method="unequal", k1=k1.origin, k2=k2.origin)
  return(re)
}



##########################################################################
#########         EMPIRICAL METHOD
##########################################################################
#########      3) COST-MINIMISING THRESHOLD (EMPIRICAL) FUNCTION
##########################################################################
##### arguments: k1=vector containing the healthy sample values 
#####		   	     k2=vector containing the diseased sample values
#####  	   	     rho=disease prevalence
#####		   	     costs=cost matrix
#####            extra.info=should extra information be added to the output? Default, FALSE
##### value: the threshold estimated by the two samples through the empirical estimator
##########################################################################
thresEmp2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), extra.info=FALSE){
  k1.origin <- k1
  k2.origin <- k2
  rho.origin <- rho
  costs.origin <- costs
  if(mean(k1)>mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  # cost matrix
  c.t.pos <- costs[1,1]
  c.t.neg <- costs[1,2]
  c.f.pos <- costs[2,1]
  c.f.neg <- costs[2,2]
  
  # empirical estimation
  n1 <- length(k1)
  n2 <- length(k2)
  n <- n1+n2
  v <- c(k1,k2)
  ind.origin <- c(rep(0,n1), rep(1,n2))
  vi <- cbind(v,ind.origin)
  ord.v <- vi[order(vi[, 1], vi[, 2]), 1:2]
  sens <- rep(NA, n)
  spec.c <- rep(NA, n) 
  # sensitivity and 1-specificity for each observation
  for(j in 1:n){
    sens[j] <- (sum(ord.v[(j:n), 2]))/n2 # Sens
    spec.c[j] <- ((n + 1 - j) - sum(ord.v[(j:n), 2]))/n1 # 1 - Spec
  }
  # cost for each observation
  cost.non.par <- n*(c.t.pos*sens*rho + c.f.neg*(1-sens)*rho + c.f.pos*spec.c*(1-rho) + c.t.neg*(1-spec.c)*(1-rho))
  # together
  total <- data.frame(ord.v, cost.non.par, sens, 1-spec.c)
  colnames(total)[5] <- "spec"
  # search and remove repeated observations
  which.duplicated <- which(duplicated(total[, "v"]))
  if (length(which.duplicated) > 0){
    total <- total[-which.duplicated, ]
  } 
  # minimum cost
  ind.min.cost <- which(total[, "cost.non.par"]==min(total[, "cost.non.par"]))
  sens.min <- total[ind.min.cost, "sens"]
  spec.min <- total[ind.min.cost, "spec"]
  cut.min <- total[ind.min.cost, "v"]
  cost.min <- total[ind.min.cost, "cost.non.par"]
  # if there are two observations that lead to the same cost...
  howmany <- length(cut.min)
  if (howmany>1){
    interval <- subset(total, v >= cut.min[1] & v<= cut.min[length(cut.min)])
    cut.min <- mean(interval$v)
    # sens, spec
    sens.min <- (sum(vi[, "v"]>cut.min & vi[, "ind.origin"]==1))/n2
    spec.min <- (sum(vi[, "v"]<cut.min & vi[, "ind.origin"]==0))/n1
    # cost
    cost.min <- n*(c.t.pos*sens.min*rho + c.f.neg*(1-sens.min)*rho + c.f.pos*(1-spec.min)*(1-rho) + c.t.neg*(spec.min)*(1-rho))
    warning(paste(howmany, "observations lead to the minimum cost. The mean of the values between them is returned as threshold."), call.=FALSE)
  }
  
  # slope
  beta <- slope(rho, costs)
  # results
  re <- list(thres=cut.min, sens=sens.min, spec=spec.min, cost=cost.min, costs=costs.origin, R=beta, prev=rho.origin, method="empirical", k1=k1.origin, k2=k2.origin)
  # for plotCostROC if extra.info is TRUE
  if (extra.info){
    re$tot.thres <- total[, "v"]
    re$tot.cost <- total[, "cost.non.par"]
    re$tot.spec.c <- spec.c
    re$tot.sens <- sens
  }
  return(re)
}

#############################################################
## Print function for class "thres2"
#############################################################
print.thres2 <- function(x, ...){
  if (x$T$method == "parametric"){
    cat("\nEstimate:")
    cat("\n  Threshold: ", x$T$thres)
    cat("\n")
    cat("\nConfidence intervals (parametric bootstrap):")
    cat("\n  CI based on normal distribution:", x$CI$low.norm, " - ", x$CI$up.norm)
    cat("\n  CI based on percentiles:", x$CI$low.perc, " - ", x$CI$up.perc)
    cat("\n  Bootstrap resamples:", x$CI$B)
    cat("\n")
    cat("\nParameters used:")
    cat("\n  Disease prevalence:", x$T$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$T$costs)
    cat("\n  R:", x$T$R)
    cat("\n  Significance Level: ", x$CI$alpha)
    cat("\n  Method:", x$T$method)
    cat("\n  Distribution assumed for the healthy sample: ", x$T$dist1, "(", round(x$T$pars1[1], 2), ", ", round(x$T$pars1[2], 2), ")", sep="")
    cat("\n  Distribution assumed for the diseased sample: ", x$T$dist2, "(", round(x$T$pars2[1], 2), ", ", round(x$T$pars2[2], 2), ")", sep="")
    cat("\n")
  }
  if (x$T$method == "empirical"){
    cat("\nEstimate:")
    cat("\n  Threshold: ", x$T$thres)
    cat("\n  Minimum Cost: ", x$T$cost)
    cat("\n")
    cat("\nConfidence intervals (bootstrap):")
    cat("\n  CI based on normal distribution:", x$CI$low.norm, " - ", x$CI$up.norm)
    cat("\n  CI based on percentiles:", x$CI$low.perc, " - ", x$CI$up.perc)
    cat("\n  Bootstrap resamples:", x$CI$B)
    cat("\n")
    cat("\nParameters used:")
    cat("\n  Disease prevalence:", x$T$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$T$costs)
    cat("\n  R:", x$T$R)
    cat("\n  Method:", x$T$method)
    cat("\n  Significance Level:", x$CI$alpha)
    cat("\n")
  }
  if (x$T$method == "equal" | x$T$method == "unequal"){
    cat("\nEstimate:")
    cat("\n  Threshold: ", x$T$thres)
    cat("\n")
    if(x$CI$ci.method == "delta"){
      cat("\nConfidence interval (delta method):")
      cat("\n  Lower Limit:", x$CI$lower)
      cat("\n  Upper Limit:", x$CI$upper)
      cat("\n")
    }
    if(x$CI$ci.method == "boot"){
      cat("\nConfidence intervals (bootstrap):")
      cat("\n  CI based on normal distribution: ", x$CI$low.norm, " - ", x$CI$up.norm)
      cat("\n  CI based on percentiles: ", x$CI$low.perc, " - ", x$CI$up.perc)
      cat("\n  Bootstrap resamples: ", x$CI$B)
      cat("\n")
    }    
    cat("\nParameters used:")
    cat("\n  Disease prevalence:", x$T$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$T$costs)
    cat("\n  R:", x$T$R)
    cat("\n  Method:", x$T$method)
    cat("\n  Significance Level: ", x$CI$alpha)
    cat("\n")
  }
}


##########################################################################
#########         GETTING PARAMETERS OF A DISTRIBUTION THROUGH fitdistr() [library MASS]
##########################################################################
getParams <- function(k, dist){
  if (dist %in% c("cauchy", "gamma", "weibull")){
    pars <- fitdistr(k, dist)$estimate
  }else if(dist=="beta"){
    # needs initial values
    sigma2 <- var(k)
    mu <- mean(k)
    shape1.start <- ((1-mu)/sigma2-1/mu)*mu^2
    shape2.start <- shape1.start*(1/mu-1)
    pars <- fitdistr(k, "beta", start=list(shape1=shape1.start, shape2=shape2.start))$estimate
  }else if(dist=="chisq"){
    # needs initial values
    sigma2 <- var(k)
    mu <- mean(k)
    ncp.start <- sigma2/2-mu
    df.start <- mu-ncp.start
    pars <- fitdistr(k, "chi-squared", start=list(df=df.start, ncp=ncp.start))$estimate
  }else if(dist=="lnorm"){
    pars <- fitdistr(k, "lognormal")$estimate
  }else if(dist=="logis"){
    pars <- fitdistr(k, "logistic")$estimate
  }else if(dist=="norm"){
    pars <- fitdistr(k, "normal")$estimate
  }
  out <- pars
  return(out)
}


##########################################################################
#########         THRESHOLD COMPUTATION
##########################################################################
##### arguments: k1=vector containing the healthy sample values 
#####  	   	     k2=vector containing the diseased sample values
#####  	   	     rho=disease prevalence
#####		   	     costs=cost matrix
#####            method=method used in the estimation. The user can choose between:
#####              "equal": assumes binormality and equal variances
#####              "unequal": assumes binormality and unequal variances
#####              "empirical": leaves out any distributional assumption
#####              "parametric": based on the distribution assumed for the two populations.
#####              The user must specify these distributions through parameters 'dist1'
#####              (for the healthy population) and 'dist2' (for the healthy population).
#####              They can be chosen between the following 2-parameter distributions:
#####              "beta", "cauchy", "chisq" (chi-squared), "gamma", "lnorm" (lognormal),
#####              "logis" (logistic), "norm" (normal) and "weibull".
#####              Its parameters are estimated from the samples k1 and k2. Uses 'thresTH2()'.
#####            dist1, dist2: distribution to be assumed for healthy and disease populations,
#####            respectively. Ignored when method!="parametric".
#####            Default, "equal". The user can specify just the initial letters.
#####            ci.method=method to be used for the confidence intervals computation. The user can
#####            choose between:
#####              "delta": delta method is used to estimate the threshold standard error
#####              assuming a binormal underlying model.
#####              "boot": the confidence interval is computed by bootstrap.
#####            Default, "delta". The user can specify just just the initial letters.
#####            B=number of bootstrap resamples when ci.method="boot". Otherwise, ignored.
#####            Default, 1000.
#####  	   	     alpha=significance level for the confidence interval. Default, 0.05.
#####            extra.info=when using method="empirical", if set to T it returns
#####            extra information about the computation of the threshold. Ignored when method!="empirical".
#####            na.rm=a logical value indicating whether NA values in k1 and k2 should
#####            be stripped before the computation proceeds. Default, FALSE.
##### value: the threshold estimated
##########################################################################
thres2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), method=c("equal", "unequal", "empirical", "parametric"), dist1=NULL, dist2=NULL, ci.method=c("delta", "boot"), B=1000, alpha=0.05, extra.info=FALSE, na.rm=FALSE){
  # error handling
  if (!(rho > 0 & rho < 1)){
    stop("The disease prevalence rho must be a number in (0,1)")
  }
  if (!is.matrix(costs)){
    stop("'costs' must be a matrix")
  }
  if (dim(costs)[1] != 2 | dim(costs)[2] != 2){
    stop("'costs' must be a 2x2 matrix")
  }
  if (!is.numeric(k1) | !is.numeric(k2)){
    stop("'k1' and 'k2' must be numeric vectors")
  }
  # NAs handling
  if (na.rm){
    k1 <- k1[!is.na(k1)]
    k2 <- k2[!is.na(k2)]
  }
  # function
  method <- match.arg(method)
  ci.method <- match.arg(ci.method)
  if (method=="equal"){
    T <- thresEq2(k1, k2, rho, costs)
    if (ci.method=="delta"){
      ci <- icDeltaEq2(k1, k2, rho, costs, T$thres, a=alpha)
    }
    if (ci.method=="boot"){
      ci <- icBootEq2(k1, k2, rho, costs, T$thres, B=B, a=alpha)
    }
  }
  if (method=="unequal"){
    T <- thresUn2(k1, k2, rho, costs)
    if (ci.method=="delta"){
      ci <- icDeltaUn2(k1, k2, rho, costs, T$thres, a=alpha)
    }
    if (ci.method=="boot"){
      ci <- icBootUn2(k1, k2, rho, costs, T$thres, B=B, a=alpha)
    }
  }
  if (method=="empirical"){
    if (ci.method=="delta"){
      stop("When method='empirical', CIs cannot be computed based on delta method (choose ci.method='boot')")
    }
    T <- thresEmp2(k1, k2, rho, costs, extra.info)
    ci <- icEmp2(k1, k2, rho, costs, T$thres, B=B, a=alpha)
  }
  if (method=="parametric"){
    # error handling
    if (ci.method=="delta"){
      stop("When method='parametric', CIs cannot be computed based on delta method (choose ci.method='boot')")
    }
    if (is.null(dist1) | is.null(dist2)){
      stop("When method='parametric', 'dist1' and 'dist2' must be specified")
    }
    if (dist1=="norm" & dist2=="norm"){
      stop("When assuming a binormal distribution, choose method='equal' or 'unequal'")
    }
    if (!(dist1 %in% c("beta", "cauchy", "chisq", "gamma", "lnorm", "logis", "nbinom", "norm", "weibull"))){
      stop("Unsupported distribution for 'dist1'")
    }
    if (!(dist2 %in% c("beta", "cauchy", "chisq", "gamma", "lnorm", "logis", "nbinom", "norm", "weibull"))){
      stop("Unsupported distribution for 'dist2'")
    }
    # parameter estimation through 'fitdistr()'
    # dist1
    pars1 <- getParams(k1, dist1)
    # dist2
    pars2 <- getParams(k2, dist2)
    # threshold+CI estimation
    T <- thresTH2(dist1, dist2, pars1[1], pars1[2], pars2[1], pars2[2], rho, costs)
    T <- unclass(T)
    T$method <- "parametric"
    # adding some information to the output
    T$k1 <- k1
    T$k2 <- k2
    T$dist1 <- dist1
    T$dist2 <- dist2
    T$pars1 <- pars1
    T$pars2 <- pars2
    # confidence interval by parametric bootstrap
    ci <- icBootTH(dist1, dist2, pars1[1], pars1[2], pars2[1], pars2[2], length(k1), length(k2), rho, costs, T$thres, B=B, a=alpha)
  }
  out <- list(T=T, CI=ci)
  class(out) <- "thres2"
  return(out)
}


##########################################################################
############      VARIANCE ESTIMATORS
############	    EQUAL VARIANCES
##########################################################################
############       FUNCTION OF ML VARIANCE OF VARIANCE ESTIMATOR
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
##### value: Maximum-Likelihood variance of variance estimator (equal variances assumed)
##########################################################################
varVarEq <- function(k1, k2){
  n1 <- length(k1)
  n2 <- length(k2)
  est <- 2*(varPooled(k1,k2))^2/(n1 + n2 - 1)
  return(est)
}
##########################################################################
############       FUNCTION OF Maximum Likelihood VARIANCE OF MEAN ESTIMATOR  (DIS+NON.DIS)
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      t=size of the sample whose mean estimator variance we wish to estimate 
##### value: Maximum-Likelihood variance of mean estimator (equal variances assumed)
##########################################################################
varMeanEq <- function(k1, k2, t){
  est <- varPooled(k1, k2)/t
  return(est)
}

##########################################################################
##########################################################################
############       PARTIAL DERIVATIVES... (Skaltsa et al 2010)
##########################################################################
############       ...OF COMMON VARIANCE
##### arguments:  k1=vector containing the healthy sample values 
#####		          k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: threshold estimator partial derivative with respect to the diseased or healthy population
#####              variance
##########################################################################
derVarEq <- function(k1, k2, rho, costs){
  est <- log(slope(rho,costs))/(mean(k2) - mean(k1))
  return(est)
}

##########################################################################
############       ...OF DISEASED MEAN
##### arguments:  k1=vector containing the healthy sample values 
#####		        	k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: threshold estimator partial derivative with respect to the diseased population mean
##########################################################################
derMeanDisEq <- function(k1, k2, rho, costs){
  est <- 1/2 - (varPooled(k1,k2)*log(slope(rho,costs)))/((mean(k2) - mean(k1))^2)
  return(est)
}

###################################################################################
############      ...OF NON-DISEASED MEAN
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: threshold estimator partial derivative with respect to the diseased population variance
##################################################################################
derMeanNDisEq <- function(k1, k2, rho, costs){
  est <- 1/2 + (varPooled(k1,k2)*log(slope(rho,costs)))/((mean(k2) - mean(k1))^2)
  return(est)
}

##################################################################################
############ 1) FUNCTION OF VARIANCE OF THRESHOLD (EQUAL VARIANCES)
############       DELTA METHOD
##### arguments:  k1=vector containing the healthy sample values 
#####		        	k2=vector containing the diseased sample values 
#####		         	rho=diseased sample size
#####		  	      costs=cost matrix
##### value: delta method threshold variance (assuming equal variances)
##################################################################################
varDeltaEq2 <- function(k1, k2, rho, costs){
  if(mean(k1) > mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  # n's
  n1 <- length(k1)
  n2 <- length(k2)
  # matrices
  d <- matrix(c(derMeanNDisEq(k1, k2, rho, costs), derMeanDisEq(k1, k2, rho, costs), derVarEq(k1, k2, rho, costs)), byrow=T, ncol=3, nrow=1)
  sigma <- matrix(c(varMeanEq(k1, k2, n2), 0, 0, 0, varMeanEq(k1, k2, n1), 0, 0, 0, varVarEq(k1, k2)), byrow=T, ncol=3, nrow=3)
  # Var(T)=d*sigma*d^T, Skaltsa et al 2010
  est <- d%*%sigma%*%t(d)
  return(est)
}

##################################################################################
############	 CONFIDENCE INTERVAL FUNCTION 1 - EQUAL VARIANCES
############       DELTA METHOD-NORMAL APPROXIMATION
##### arguments:  	k1=vector containing the healthy sample values 
#####		   	        k2=vector containing the diseased sample values 
#####		   	        rho=diseased sample size
#####		  	        costs=cost matrix
#####               Thres=threshold point estimate
#####		   	        a=significance level
##### value: confidence interval limits
##################################################################################
icDeltaEq2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), Thres, a=0.05){
  stdev <- sqrt(varDeltaEq2(k1, k2, rho, costs))  
  ic1 <- Thres + qnorm(a/2)*stdev
  ic2 <- Thres + qnorm(1-a/2)*stdev  
  ic <- list(lower=ic1, upper=ic2, alpha=a, ci.method="delta")
  return(ic)
}





##################################################################################
############       VARIANCES OF ESTIMATORS
############       UNEQUAL VARIANCES
##################################################################################
############       FUNCTION OF ML VARIANCE OF VARIANCE ESTIMATOR
##### arguments:  k=vector containing the sample values
##### 		        t=sample size
##### value: Maximum-Likelihood variance of variance estimator (unequal variances assumed)
##################################################################################
varVarUn <- function(k, t){
  est <- 2*(var(k))^2/(t - 1)
  return(est)
} 
##################################################################################
############       FUNCTION OF ML VARIANCE OF MEAN ESTIMATOR  (DIS+NON.DIS)
##### arguments:  k=vector containing the sample values
##### 		        t=sample size
##### value: Maximum-Likelihood variance of mean estimator (unequal variances assumed)
##################################################################################
varMeanUn <- function(k, t){
  est <- var(k)/t
  return(est)
}

#################################################################################
############       PARTIAL DERIVATIVES...
#################################################################################
############       ...OF DISEASED MEAN
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      rho=diseased sample size
#####		  	      costs=cost matrix
##### value: threshold estimator partial derivative with respect to the diseased population mean
#####              (assuming unequal variances)
#################################################################################
derMeanDisUn <- function(k1, k2, rho, costs){
  est <- (sd(k2)*sd(k1)*(mean(k2) - mean(k1))/sqrt((mean(k2) - mean(k1))^2+2*(var(k2) - var(k1))*log(sd(k2)*slope(rho,costs)/sd(k1))) - var(k1))/(var(k2) - var(k1)) 
  return(est)
}

#################################################################################
############       ...OF NON-DISEASED MEAN
##### arguments:  k1=vector containing the healthy sample values 
#####		   	k2= vector containing the diseased sample values 
#####		   	rho=diseased sample size
#####		  	costs=cost matrix
##### value: threshold estimator partial derivative with respect to the healthy population mean
#####              (assuming unequal variances)
#################################################################################
derMeanNDisUn <- function(k1, k2, rho, costs){
  est <- (var(k2) - sd(k2)*sd(k1)*(mean(k2) - mean(k1))/sqrt((mean(k2) - mean(k1))^2+2*(var(k2) - var(k1))*log(sd(k2)*slope(rho,costs)/sd(k1))))/(var(k2) - var(k1))
  return(est)
}

#################################################################################
############       ...OF DISEASED VARIANCE
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: threshold estimator partial derivative with respect to the diseased population variance
#####              (assuming unequal variances)
#################################################################################
derVarDisUn <- function(k1, k2, rho, costs){
  beta <- slope(rho,costs)
  est <- (sd(k1)*sd(k2)*((var(k2) - var(k1))/var(k2) + 2*log(beta*sd(k2)/sd(k1)))/(2*sqrt(2*log(beta*sd(k2)/sd(k1))*(var(k2) - var(k1))+(mean(k2) - mean(k1))^2)) + sd(k1)*sqrt(2*log(beta*sd(k2)/sd(k1))*(var(k2) - var(k1)) + (mean(k2) - mean(k1))^2)/(2*sd(k2)) + mean(k1)) /(var(k2) - var(k1)) - (sd(k1)*sd(k2)*sqrt(2*log(beta*sd(k2)/sd(k1))*(var(k2) - var(k1))+(mean(k2)-mean(k1))^2)+mean(k1)*var(k2)-mean(k2)*var(k1))/(var(k2)-var(k1))^2 
  return(est)
}

#################################################################################
############       ...OF NON-DISEASED VARIANCE
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: threshold estimator partial derivative with respect to the healthy population variance
#####              (assuming unequal variances)
#################################################################################
derVarNDisUn <- function(k1, k2, rho, costs){
  beta <- slope(rho, costs)
  est <- (-mean(k2)*var(k1)+sd(k2)*sqrt(2*log(sd(k2)*beta/sd(k1))*(var(k2)-var(k1))+(mean(k2)-mean(k1))^2)*sd(k1)+mean(k1)*var(k2)) /(var(k2)-var(k1))^2 +(sd(k2)*(-(var(k2)-var(k1))/var(k1)-2*log(sd(k2)*beta/sd(k1)))*sd(k1) /(2*sqrt(2*log(sd(k2)*beta/sd(k1))*(var(k2)-var(k1))+(mean(k2)-mean(k1))^2)) +sd(k2)*sqrt(2*log(sd(k2)*beta/sd(k1))*(var(k2)-var(k1))+(mean(k2)-mean(k1))^2)/(2*sd(k1))-mean(k2)) /(var(k2)-var(k1)) 
  return(est)
}

#################################################################################
############ 2) FUNCTION OF VARIANCE OF THRESHOLD (UNEQUAL VARIANCES)
############       DELTA METHOD
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: delta method threshold variance (assuming unequal variances)
#################################################################################
varDeltaUn2 <- function(k1, k2, rho, costs){
  if(mean(k1) > mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  # n's
  n1 <- length(k1)
  n2 <- length(k2)
  ctrl <- sqroot(k1,k2,rho,costs) 	
  if(ctrl < 0){
    est <- NA
    warning("Negative discriminant; cannot solve the second-degree equation")
  }else{
    # matrices
    d <- matrix(c(derMeanNDisUn(k1, k2, rho, costs), derMeanDisUn(k1, k2, rho, costs), derVarNDisUn(k1, k2, rho, costs), derVarDisUn(k1, k2, rho, costs)), byrow=T, ncol=4, nrow=1)
    sigma <- matrix(c(varMeanUn(k1, n1), 0, 0, 0, 0, varMeanUn(k2, n2), 0, 0, 0, 0, varVarUn(k1, n1), 0, 0, 0, 0, varVarUn(k2, n2)), byrow=T, ncol=4, nrow=4)
    # Var(T) following Skaltsa et al 2010
    est <- d%*%sigma%*%t(d)
  }
  return(est)
}

#################################################################################
############ 	 FUNCTION OF CONFIDENCE INTERVAL 1 - UNEQUAL
############       DELTA METHOD-NORMAL APPROXIMATION
##### arguments:  k1=vector containing the healthy sample values 
#####		         	k2=vector containing the diseased sample values 
#####		   	      rho=diseased sample size
#####		  	      costs=cost matrix
#####		   	      a=significance level
##### value: confidence interval limits
##################################################################################
icDeltaUn2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), Thres, a=0.05){
  stdev <- sqrt(varDeltaUn2(k1, k2, rho, costs))	
  ic1 <- Thres + qnorm(a/2)*stdev
  ic2 <- Thres + qnorm(1-a/2)*stdev
  ic <- list(lower=ic1, upper=ic2, alpha=a, ci.method="delta")  
  return(ic)
}





################################################################################
############       DATA GENERATION
############ 	 NON-PARAMETRIC RESAMPLING
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      B=number of bootstrap resamples
##### value: two-object list [[1]]:healthy resample matrix, [[2]]:diseased resample matrix
################################################################################
resample2 <- function(k1, k2, B){
  n1 <- length(k1)
  n2 <- length(k2)
  t0 <- matrix(sample(k1, n1*B, replace=TRUE), nrow=n1) # B resamples of k1 with replacement
  t1 <- matrix(sample(k2, n2*B, replace=TRUE), nrow=n2) # B resamples of k2 with replacement
  t <- list(t0, t1)
  return(t)
}

###################################################################################
###################################################################################
############          CONFIDENCE INTERVALS
###################################################################################
############               BOOTSTRAP
###################################################################################
###################################################################################
###############    FUNCTION OF BOOTSTRAP CONFIDENCE INTERVALS 
###############    EQUAL VARIANCES
##########              1) BOOTSTRAP VERSION OF S.E.- NORMAL
##########              2) PERCENTILE BIAS-CORRECTED
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      B=number of bootstrap resamples
#####		   	      a=significance level
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### returns: the 2 interval limits
##################################################################################
icBootEq2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), Thres, B, a=0.05){
  if(mean(k1) > mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  # resamples
  t <- resample2(k1,k2,B)
  t0 <- t[[1]]
  t1 <- t[[2]]
  
  # bootstrap
  cut <- sapply(1:B,function(i){thresEq2(t0[,i],t1[,i],rho,costs)[[1]]})
  est.se <- sd(cut)
  
  ###### 1) NORMAL-BOOTSTRAP SE
  norm.bootSE <- c(Thres + qnorm(a/2)*est.se, Thres + qnorm(1-a/2)*est.se)
  
  ###### 2) PERCENTILE
  percentil <- (c(quantile(cut,a/2), quantile(cut,1-a/2)))
  
  # results
  re <- list(low.norm=norm.bootSE[1], up.norm=norm.bootSE[2], low.perc=percentil[1], up.perc=percentil[2], alpha=a, B=B, ci.method="boot")
  return(re)
}


###################################################################################
###################################################################################
###############     BOOTSTRAP CONFIDENCE INTERVAL FUNCTION
###############     UNEQUAL VARIANCES
##########               1) BOOTSTRAP VERSION OF S.E.- NORMAL
##########               2) PERCENTILE
##### arguments:   k1=vector containing the healthy sample values 
#####		   	       k2=vector containing the diseased sample values 
#####		   	       B=number of bootstrap resamples
#####		   	       a=significance level
#####		   	       rho=disease prevalence
#####		  	       costs=cost matrix
##### value: the 2 interval limits
##################################################################################
icBootUn2 <- function(k1, k2, rho, costs=matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE), Thres, B, a=0.05){
  if(mean(k1) > mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  
  t <- resample2(k1,k2,B)
  t0 <- t[[1]]
  t1 <- t[[2]]
  
  cut <- sapply(1:B,function(i){thresUn2(t0[,i],t1[,i],rho,costs)[[1]]})
  
  est.se <- sd(na.omit(cut))
  
  ###### 1) NORMAL-BOOTSTRAP SE
  norm.bootSE <- c(Thres + qnorm(a/2)*est.se, Thres + qnorm(1-a/2)*est.se)
  
  ###### 2) PERCENTILE
  percentil <- (c(quantile(na.omit(cut),a/2), quantile(na.omit(cut),1-a/2)))
  
  # results
  re <- list(low.norm=norm.bootSE[1], up.norm=norm.bootSE[2], low.perc=percentil[1], up.perc=percentil[2], alpha=a, B=B, ci.method="boot")
  return(re)
}





##################################################################################
##################################################################################
################    BOOTSTRAP CONFIDENCE INTERVAL
################    EMPIRICAL
##################################################################################
############           1) BOOTSTRAP VERSION OF S.E.- NORMAL
############           2) BOOTSTRAP PERCENTILE
##### arguments:  k1=vector containing the healthy sample values 
#####		   	      k2=vector containing the diseased sample values 
#####		   	      B=number of bootstrap resamples
#####		   	      a=significance level
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### value: the 2 interval limits
##################################################################################
icEmp2 <- function(k1,k2,rho,costs=matrix(c(0,0,1,(1-rho)/rho),2,2, byrow=TRUE),Thres,B=500,a=0.05){
  if(mean(k1)>mean(k2)){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- k1; k1 <- k2; k2 <- g
  }
  
  t <- resample2(k1,k2,B)
  t0 <- t[[1]]
  t1 <- t[[2]]
  
  cut <- rep(NA,B)	
  cut <- suppressWarnings(sapply(1:B,function(j){thresEmp2(t0[,j],t1[,j],rho,costs)[[1]]}))
  est.se <- sd(cut)
  
  ###### 1) NORMAL-BOOTSTRAP SE
  norm <- c(Thres+qnorm(a/2)*est.se, Thres+qnorm(1-a/2)*est.se)
  
  ###### 2) PERCENTILE
  percentil <- (c(quantile(cut,a/2), quantile(cut,1-a/2)))
  
  # results
  re <- list(low.norm=norm[1], up.norm=norm[2], low.perc=percentil[1], up.perc=percentil[2], alpha=a, B=B, ci.method="boot")
  return(re)
}

################################################################################
############       DATA GENERATION
############    PARAMETRIC RESAMPLING
##### arguments:  dist1, dist2=distributions assumed for the healthy and diseased population, respectively. 
#####  	   	      B=number of bootstrap resamples
##### value: two-object list [[1]]:healthy resample matrix, [[2]]:diseased resample matrix
################################################################################
aux.par.boot <- function(dist1, dist2, par1.1, par1.2, par2.1, par2.2, n1, n2, B){
  t0 <- matrix(rand(dist1)(n1*B, par1.1, par1.2), nrow=n1) # B resamples of dist1(par1.1, par1.2) with replacement
  t1 <- matrix(rand(dist2)(n2*B, par2.1, par2.2), nrow=n2) # B resamples of dist2(par2.1, par2.2) with replacement
  t <- list(t0,t1)
  return(t)
}

###################################################################################
###################################################################################
###############    FUNCTION OF PARAMETRIC BOOTSTRAP FOR THEORETICAL METHOD
##########              1) BOOTSTRAP VERSION OF S.E.- NORMAL
##########              2) PERCENTILE BIAS-CORRECTED
##### arguments:  dist1, dist2=distributions assumed for the healthy and diseased population, respectively. 
#####             par1.1, par1.2, par2.1, par2.2, n1, n2=parameters and sample sizes
#####		   	      B=number of bootstrap resamples
#####		   	      a=significance level
#####		   	      rho=disease prevalence
#####		  	      costs=cost matrix
##### returns: the 2 interval limits
##################################################################################
icBootTH <- function(dist1, dist2, par1.1, par1.2, par2.1, par2.2, n1, n2, rho, costs=matrix(c(0,0,1,(1-rho)/rho), 2, 2, byrow=TRUE), Thres, B=500, a=0.05){
  median1 <- quant(dist1)(0.5, par1.1, par1.2)
  median2 <- quant(dist2)(0.5, par2.1, par2.2)
  if(median1 > median2){
    rho <- 1-rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- par2.1; par2.1 <- par1.1; par1.1 <- g
    f <- par2.2; par2.2 <- par1.2; par1.2 <- f
    auxdist <- dist2; dist2 <- dist1; dist1 <- auxdist
  }
  # bootstrap resamples
  t <- aux.par.boot(dist1, dist2, par1.1, par1.2, par2.1, par2.2, n1, n2, B)
  t0 <- t[[1]]
  t1 <- t[[2]]
  # computing threshold of the resamples...
  # 1) parameter estimation for each resample
  pars1 <- sapply(1:B, function(i){getParams(t0[, i], dist1)})
  pars2 <- sapply(1:B, function(i){getParams(t1[, i], dist2)})
  # 2) TH cut
  cut <- sapply(1:B,function(i){thresTH2(dist1, dist2, pars1[1, i], pars1[2, i], pars2[1, i], pars2[2, i], rho, costs)[[1]]})
  # sd
  est.se <- sd(na.omit(cut))
  
  ###### 1) NORMAL-BOOTSTRAP SE
  norm.bootSE <- c(Thres + qnorm(a/2)*est.se, Thres + qnorm(1-a/2)*est.se)
  
  ###### 2) PERCENTILE
  percentil <- (c(quantile(na.omit(cut),a/2), quantile(na.omit(cut),1-a/2)))
  
  # results
  re <- list(low.norm=norm.bootSE[1], up.norm=norm.bootSE[2], low.perc=percentil[1], up.perc=percentil[2], alpha=a, B=B, ci.method="boot")
  return(re)
}



##################################################################################
###### 	PLOTCOSTROC ("thres2" and "thres3" classes)
##################################################################################
######		COST FUNCTION PLOT + ROC CURVE, OPTIMAL THRESHOLD, SENSITIVITY, SPECIFICITY (2-state)
######  2-states: plotCost function provides two graphics:
######  1) the cost function with the cost minimising threshold in red
######  2) the ROC curve with the sensitivity and specificity achieved in red
######  3-states: plotCost function provides two graphics:
######  The contribution of both thresholds to the cost function, with the minimising thresholds in res
######  arguments:
######        x: 'thres2' or 'thres3' object 
######        type: type of pline for plots
##################################################################################
plotCostROC <- function(x, type="l", ...){
  if (!(class(x) %in% c("thres2", "thres3"))){
    stop("'x' must be a 'thres2' or 'thres3' object")
  }
  if (class(x)=="thres2"){ # 2 states
    if (x$T$method == "empirical"){
      if (length(x$T) != 14){
        stop("use argument 'extra.info = TRUE' in 'thres2()'")
      }
      # plots
      # plot 1: empirical cost function
      plot(x$T$tot.thres, x$T$tot.cost, xlab="Threshold", ylab="Cost", main="Empirical Cost function", type=type, ...)
      points(x$T$thres, x$T$cost, col=2, pch=19)
      par(ask=T)
      # plot 2: empirical ROC curve
      plot(x$T$tot.spec.c,x$T$tot.sens, main="Empirical ROC curve", xlab="1-Specificity", ylab="Sensitivity", type=type, ...)
      points(1-x$T$spec, x$T$sens, col=2, pch=19)
      par(ask=F)
    }
    if (x$T$method =="equal"){
      k1 <- x$T$k1; k2 <- x$T$k2; rho <- x$T$prev; costs <- x$T$costs
      # changes
      changed <- FALSE
      if (mean(k1)>mean(k2)){
        rho <- 1-rho
        costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
        g <- k1; k1 <- k2; k2 <- g
        changed <- TRUE
      }
      # plot limits
      k <- c(x$T$k1, x$T$k2)
      rangex <- range(k)
      xlim <- c(floor(rangex[1]), ceiling(rangex[2]))
      # cost function
      cost <- function(t){
        sd <- sqrt(varPooled(k1, k2))
        TP <- (1-pnorm(t, mean(k2), sd))*rho
        FN <- (pnorm(t, mean(k2), sd))*rho
        FP <- (1-pnorm(t, mean(k1), sd))*(1-rho)
        TN <- (pnorm(t, mean(k1), sd))*(1-rho)
        n <- length(k1)+length(k2)
        C <- n*(TP*costs[1,1]+FN*costs[2, 2]+FP*costs[2, 1]+TN*costs[1, 2])
        return(as.numeric(C))
      }
      # plot 1: cost function
      plot(cost, xlim=xlim, type=type, xlab="t", ylab="cost(t)", ...)
      points(x$T$thres, cost(x$T$thres), col=2, pch=19)
      par(ask=T)
      # plot 2: ROC curve
      CUT <- x$T$thres
      if (!changed){
        resp.CUT <- ifelse(k<CUT, 0, 1)
      }else{
        resp.CUT <- ifelse(k>CUT, 0, 1)
      }   
      resp <- c(rep(0, length(x$T$k1)), rep(1, length(x$T$k2)))
      resp.CUT <- factor(resp.CUT, c("0", "1"))
      roc <- roc(response=resp, predictor=k)
      plot(roc)
      # add sens and spec given by the threshold
      tab <- table(resp.CUT, resp)[2:1, 2:1]
      SENS <- tab[1, 1]/(tab[1, 1]+tab[2, 1])
      SPEC <- tab[2, 2]/(tab[2, 2]+tab[1, 2])
      points(SPEC, SENS, col=2, pch=19)       
      par(ask=F)    
    }
    if (x$T$method=="unequal"){
      k1 <- x$T$k1; k2 <- x$T$k2; rho <- x$T$prev; costs <- x$T$costs
      # changes
      changed <- FALSE
      if (mean(k1)>mean(k2)){
        rho <- 1-rho
        costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
        g <- k1; k1 <- k2; k2 <- g
        changed <- TRUE
      }
      # plot limits
      k <- c(x$T$k1, x$T$k2)
      rangex <- range(k)
      xlim <- c(floor(rangex[1]), ceiling(rangex[2]))
      # cost function
      cost <- function(t){
        sd1 <- sd(k1)
        sd2 <- sd(k2)
        TP <- (1-pnorm(t, mean(k2), sd2))*rho
        FN <- (pnorm(t, mean(k2), sd2))*rho
        FP <- (1-pnorm(t, mean(k1), sd1))*(1-rho)
        TN <- (pnorm(t, mean(k1), sd1))*(1-rho)
        n <- length(k1)+length(k2)
        C <- n*(TP*costs[1,1]+FN*costs[2, 2]+FP*costs[2, 1]+TN*costs[1, 2])
        return(as.numeric(C))
      }
      # plot 1: cost function
      plot(cost, xlim=xlim, type=type, xlab="t", ylab="cost(t)", ...)
      points(x$T$thres, cost(x$T$thres), col=2, pch=19)
      par(ask=T)
      # plot 2: ROC curve
      CUT <- x$T$thres
      if (!changed){
        resp.CUT <- ifelse(k<CUT, 0, 1)
      }else{
        resp.CUT <- ifelse(k>CUT, 0, 1)
      }   
      resp <- c(rep(0, length(x$T$k1)), rep(1, length(x$T$k2)))
      resp.CUT <- factor(resp.CUT, c("0", "1"))
      roc <- roc(response=resp, predictor=k)
      plot(roc)
      # add sens and spec given by the threshold
      tab <- table(resp.CUT, resp)[2:1, 2:1]
      SENS <- tab[1, 1]/(tab[1, 1]+tab[2, 1])
      SPEC <- tab[2, 2]/(tab[2, 2]+tab[1, 2])
      points(SPEC, SENS, col=2, pch=19)
      par(ask=F)   
    }
     if (x$T$method=="parametric"){
       k1 <- x$T$k1; k2 <- x$T$k2; rho <- x$T$prev; costs <- x$T$costs
       dist1 <- x$T$dist1; dist2 <- x$T$dist2
       par1.1 <- x$T$pars1[1]; par1.2 <- x$T$pars1[2]; par2.1 <- x$T$pars2[1]; par2.2 <- x$T$pars2[2]
       # changes
       changed <- FALSE
       if (median(k1)>median(k2)){
         rho <- 1-rho
         costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
         g <- k1; k1 <- k2; k2 <- g
         changed <- TRUE
         g <- par2.1; par2.1 <- par1.1; par1.1 <- g
         f <- par2.2; par2.2 <- par1.2; par1.2 <- f
         auxdist <- dist2; dist2 <- dist1; dist1 <- auxdist        
       }
       # plot limits
       k <- c(x$T$k1, x$T$k2)
       rangex <- range(k)
       xlim <- c(floor(rangex[1]), ceiling(rangex[2]))
       # cost function
       cost <- function(t){
         TP <- (1-p(dist2)(t, par2.1, par2.2))*rho
         FN <- (p(dist2)(t, par2.1, par2.2))*rho
         FP <- (1-p(dist1)(t, par1.1, par1.2))*(1-rho)
         TN <- (p(dist1)(t, par1.1, par1.2))*(1-rho)
         n <- length(k1)+length(k2)
         C <- n*(TP*costs[1,1]+FN*costs[2, 2]+FP*costs[2, 1]+TN*costs[1, 2])
         return(as.numeric(C))
       }
       # plot 1: cost function
       plot(cost, xlim=xlim, type=type, xlab="t", ylab="cost(t)", ...)
       points(x$T$thres, cost(x$T$thres), col=2, pch=19)
       par(ask=T)
       # plot 2: ROC curve
       CUT <- x$T$thres
       if (!changed){
         resp.CUT <- ifelse(k<CUT, 0, 1)
       }else{
         resp.CUT <- ifelse(k>CUT, 0, 1)
       }   
       resp <- c(rep(0, length(x$T$k1)), rep(1, length(x$T$k2)))
       resp.CUT <- factor(resp.CUT, c("0", "1"))
       roc <- roc(response=resp, predictor=k)
       plot(roc)
       # add sens and spec given by the threshold
       tab <- table(resp.CUT, resp)[2:1, 2:1]
       SENS <- tab[1, 1]/(tab[1, 1]+tab[2, 1])
       SPEC <- tab[2, 2]/(tab[2, 2]+tab[1, 2])
       points(SPEC, SENS, col=2, pch=19)
       par(ask=F) 
     }
  }else{ # 3-states
    k <- c(x$T$k1, x$T$k2, x$T$k3)
    rangex <- range(k)
    xlim <- c(floor(rangex[1]), ceiling(rangex[2]))
    # cost function
    dist1 <- x$T$dist1; dist2 <- x$T$dist2; dist3 <- x$T$dist3
    rho <- x$T$prev; costs <- x$T$costs
    if (dist1 =="norm" & dist2=="norm" & dist3=="norm"){
      par1.1 <- mean(x$T$k1); par1.2 <- sd(x$T$k1)
      par2.1 <- mean(x$T$k2); par2.2 <- sd(x$T$k2)
      par3.1 <- mean(x$T$k3); par3.2 <- sd(x$T$k3)
    }else{
      par1.1 <- x$T$pars1[1]; par1.2 <- x$T$pars1[2]
      par2.1 <- x$T$pars2[1]; par2.2 <- x$T$pars2[2]
      par3.1 <- x$T$pars3[1]; par3.2 <- x$T$pars3[2]
    }
    cost.t1 <- function(t){
      aux1 <- rho[1]*p(dist1)(t, par1.1, par1.2)*(costs[1, 1]-costs[1, 2])
      aux2 <- rho[2]*p(dist2)(t, par2.1, par2.2)*(costs[2, 1]-costs[2, 2])
      aux3 <- rho[3]*p(dist3)(t, par3.1, par3.2)*(costs[3, 1]-costs[3, 2])
      C <- aux1+aux2+aux3
      return(as.numeric(C))
    }
    cost.t2 <- function(t){
      aux1 <- rho[1]*p(dist1)(t, par1.1, par1.2)*(costs[1, 2]-costs[1, 3])
      aux2 <- rho[2]*p(dist2)(t, par2.1, par2.2)*(costs[2, 2]-costs[2, 3])
      aux3 <- rho[3]*p(dist3)(t, par3.1, par3.2)*(costs[3, 2]-costs[3, 3])
      C <- aux1+aux2+aux3
      return(as.numeric(C))
    }
    # plot 1: C(T1)
    plot(cost.t1, xlim=xlim, type=type, ylab="Cost(thres1)", xlab="thres1", ...)
    points(x$T$thres1, cost.t1(x$T$thres1), col=2, pch=19)
    par(ask=T)
    # plot 2: C(T2)
    plot(cost.t2, xlim=xlim, type=type, ylab="Cost(thres2)", xlab="thres2", ...)
    points(x$T$thres2, cost.t2(x$T$thres2), col=2, pch=19)
    par(ask=F)
  }
}

##################################################################################
######	DENSITY PLOT
##################################################################################
###### plot.thres2 function provides a graphic including the sample densities (diseased
###### and non-diseased populations), the threshold and its confidence interval
###### arguments:
######  x: 'thres2' object
######  bw: vector containing the bandwith for the non-diseased sample in the 1st
######      position and the bandwith for the diseased sample in the 2nd position
######      (to be passed to 'density()'). Default, c("nrd0", "nrd0").
######  ci: should the confidence interval be plotted? Default, T.
######  which.boot: in case 'x' contains CI computed by bootstrapping, which one should be printed?
######              the user can choose between "norm" (based on normal distribution)
######              or "perc" (based on percentiles). Default, "norm". This argument
######              is ignored if the CI were computed by the delta method.
######  col: 3-dimensional vector:
######       col[1]: color for the density of the non-diseased sample
######       col[2]: color for the density for the diseased sample
######       col[3]: color for the threshold and its corresponding CI
######       Default, c(1, 2, 1).
######  lty: 4-dimensional vector:
######       lty[1]: line type for the density of the non-diseased sample
######       lty[2]: line type for the density of the diseased sample
######       lty[3]: line type for the threshold
######       lty[4]: line type for the CI
######       Default, c(1, 1, 1, 2).
######  lwd: 3-dimensional vector:
######       lwd[1]: line width for the density of the non-diseased sample
######       lwd[2]: line width for the density for the diseased sample
######       lwd[3]: line width for the threshold and its corresponding CI
######       Default, c(1, 1, 1).
######  main, xlab, ...: further arguments to be passed to 'plot()'.
######  legend: logical asking if an automatic legend should be added to the graph. Default, TRUE.
######  leg.pos: position of the legend. Default, "topleft". Ignored if legend=FALSE.
######  leg.cex: a number that reescales the size of the legend. Ignored if legend=FALSE. Default, 1.
##################################################################################
plot.thres2 <- function(x, bw=c("nrd0", "nrd0"), ci=TRUE, which.boot=c("norm", "perc"), col=c(1, 2, 1), lty=c(1, 1, 1, 2), lwd=c(1, 1, 1), main=paste0("Threshold estimate ", ifelse(ci, "and CI ", ""), "(method ", x$T$method, ")"), xlab="", legend=TRUE, leg.pos="topleft", leg.cex=1, ...){
  which.boot <- match.arg(which.boot)
  k1 <- x$T$k1
  k2 <- x$T$k2
  dens.k1 <- density(k1, bw=bw[1])
  dens.k2 <- density(k2, bw=bw[2])
  min.x <- min(min(dens.k1$x), min(dens.k2$x))
  max.x <- max(max(dens.k1$x), max(dens.k2$x))
  max.y <- max(max(dens.k1$y), max(dens.k2$y))
  plot(dens.k1, xlim=c(min.x, max.x), ylim=c(0, max.y), col=col[1], lty=lty[1], lwd=lwd[1], main=main, xlab=xlab, ...)
  lines(dens.k2, col=col[2], lty=lty[2], lwd=lwd[2])
  # thres
  abline(v=x$T$thres, col=col[3], lty=lty[3], lwd=lwd[3])
  # CI
  if(ci){
    if (x$CI$ci.method != "boot"){
      abline(v=c(x$CI$lower, x$CI$upper), col=col[3], lty=lty[4], lwd=lwd[3])
    }else{
      abline(v=c(x$CI[paste0("low.", which.boot)], x$CI[paste0("up.", which.boot)]), col=col[3],  lty=lty[4], lwd=lwd[3])
    }
  }  
  # legend
  if (legend){
    legend(leg.pos, c(expression(bar(D)), "D", ifelse(ci, "Thres+CI", "Thres")), col=col, lty=lty, lwd=lwd, cex=leg.cex)
  }
}



##################################################################################
######  LINES DENSITY PLOT
##################################################################################
###### lines.thres2 function includes vertical lines for the thresholds and CIs in
###### a plot.thres2.
###### arguments:
######  x: 'thres2' object
######  ci: should the confidence interval be plotted? Default, TRUE.
######  which.boot: in case 'x' contains CI computed by bootstrapping, which one should be printed?
######              the user can choose between "norm" (based on normal distribution)
######              or "perc" (based on percentiles). Default, "norm". This argument
######              is ignored if the CI were computed by the delta method.
######  col: color for the threshold and its corresponding CI. Default, 1.
######  lty: 2-dimensional vector:
######       lty[1]: line type for the threshold
######       lty[2]: line type for the CI
######       Default, c(1, 2).
######  lwd: line width for the threshold and its corresponding CI. Default, 1.
######  ...: further arguments to be passed to 'abline()'.
##################################################################################
lines.thres2 <- function(x, ci=TRUE, which.boot=c("norm", "perc"), col=1, lty=c(1, 2), lwd=1, ...){
  which.boot <- match.arg(which.boot)
  # thres
  abline(v=x$T$thres, col=col, lty=lty[1], lwd=lwd)
  # CI
  if(ci){
    if (x$CI$ci.method != "boot"){
      abline(v=c(x$CI$lower, x$CI$upper), col=col, lty=lty[2], lwd=lwd, ...)
    }else{
      abline(v=c(x$CI[paste0("low.", which.boot)], x$CI[paste0("up.", which.boot)]), col=col,  lty=lty[2], lwd=lwd, ...)
    }
  }
}
   
                


############################################################################
####	    SAMPLE SIZE ESTIMATION
############################################################################

############################################################################
############################################################################
######### PARAMETRIC
############################################################################
######### LIMITING BETA
######### CHECKING POSITIVE SQUARE ROOT
##### arguments:
#####  par1.1=healthy population mean
#####  par1.2=healthy population standard deviation
#####  par2.1=diseased population mean
#####	 par2.2=diseased population standard deviation
#####  rho=disease prevalence
#####	 costs=cost matrix
############################################################################
control <- function(par1.1,par1.2,par2.1,par2.2,rho,costs){
  ctrl <- (par2.1-par1.1)^2+2*log((par2.2/par1.2)*slope(rho,costs))*(par2.2^2-par1.2^2)
  return(ctrl)
}

###########################################################################
############       STANDARD ERROR
###########################################################################
##########################################################################
############       VARIANCES OF ESTIMATORS
############       UNEQUAL VARIANCES
##########################################################################
############ FUNCTION OF ML VARIANCE OF VARIANCE ESTIMATOR
##########################################################################
parVarUn <- function(sdev, t){
  est <- 2*sdev^4/(t-1)
  return(est)
} 

##########################################################################
############ FUNCTION OF ML VARIANCE OF MEAN ESTIMATOR  (DIS+NON.DIS)
##########################################################################
parMeanUn <- function(sdev,t){
  est <- sdev^2/t
  return(est)
}

##########################################################################
############        PARTIAL DERIVATIVES.............
##########################################################################
############         OF DISEASED MEAN
##### arguments:
#####  par1.1=healthy population mean
#####  par1.2=healthy population standard deviation
#####  par2.1=diseased population mean
#####  par2.2=diseased population standard deviation
#####  rho=disease prevalence
#####	 costs=cost matrix
##########################################################################
parDerMeanDisUn <- function(par1.1,par1.2,par2.1,par2.2,rho,costs){
  est <- (par2.2*par1.2*(par2.1-par1.1)/sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs) )-par1.2^2)/(par2.2^2-par1.2^2) 
  return(est)
}

##########################################################################
############         OF NON-DISEASED MEAN
##### arguments:
#####  par1.1=healthy population mean
#####  par1.2=healthy population standard deviation
#####  par2.1=diseased population mean
#####  par2.2=diseased population standard deviation
#####  rho=disease prevalence
#####  costs=cost matrix
##########################################################################
parDerMeanNonDisUn <- function(par1.1,par1.2,par2.1,par2.2,rho,costs){
  est <- (par2.2^2-par2.2*par2.1*(par1.2-par1.1)/sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs)))/(par2.2^2-par2.1^2)
  return(est)
}
##########################################################################
############         OF DISEASED VARIANCE
##### arguments:
#####  par1.1=healthy population mean
#####  par1.2=healthy population standard deviation
#####  par2.1=diseased population mean
#####  par2.2=diseased population standard deviation
#####  rho=disease prevalence
#####  costs=cost matrix
##########################################################################
parDerVarDisUn <- function(par1.1,par1.2,par2.1,par2.2,rho,costs){
  est <- ((par1.2*par2.2*((par2.2^2-par1.2^2)/par2.2^2+2*log(slope(rho,costs)*par2.2/par1.2)))/(2*sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs)))+par1.2*sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs))/(2*par2.2)+par1.1)/(par2.2^2-par1.2^2)-(par1.2*par2.2*sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs))+par1.1*par2.2^2-par2.1*par1.2^2)/(par2.2^2-par1.2^2)^2 
  return(est)
}

##########################################################################
############         OF NON-DISEASED VARIANCE
##### arguments:
#####  par1.1=healthy population mean
#####  par1.2=healthy population standard deviation
#####  par2.1=diseased population mean
#####  par2.2=diseased population standard deviation
#####  rho=disease prevalence
#####  costs=cost matrix
##########################################################################
parDerVarNonDisUn <- function(par1.1,par1.2,par2.1,par2.2,rho,costs){
  est <- ((par1.2*par2.2*((par1.2^2-par2.2^2)/par2.2^2-2*log(slope(rho,costs)*par2.2/par1.2)))/(2*sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs)))+par2.2*sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs))/(2*par1.2)-par2.1)/(par2.2^2-par1.2^2)+(par1.2*par2.2*sqrt(control(par1.1,par1.2,par2.1,par2.2,rho,costs))+par1.1*par2.2^2-par2.1*par1.2^2)/(par2.2^2-par1.2^2)^2 
  return(est)
}

##########################################################################
############ SAMPLE SIZE RATIO AND MINIMUM REQUIRED SAMPLE SIZE
##### arguments:
#####  par1.1=healthy population mean
#####  par1.2=healthy population standard deviation
#####  par2.1=diseased population mean
#####  par2.2=diseased population standard deviation
#####  rho=disease prevalence
#####  costs=cost matrix
##########################################################################
SS <- function(par1.1, par1.2, par2.1, par2.2=NULL, rho, width, costs=matrix(c(0,0,1,(1-rho)/rho),2,2, byrow=TRUE), var.equal=FALSE, alpha=0.05){
  # error handling
  if (!(rho > 0 & rho < 1)){
    stop("The disease prevalence rho must be a number in (0,1)")
  }
  if (!is.matrix(costs)){
    stop("'costs' must be a matrix")
  }
  if (dim(costs)[1] != 2 | dim(costs)[2] != 2){
    stop("'costs' must be a 2x2 matrix")
  }
  if (width<=0){
    stop("'width' must be a positive number")
  }
  if (var.equal){
    if (!is.null(par2.2)){
      if (par1.2 != par2.2){
        stop("'var.equal' is set to TRUE, but 'par1.2' and 'par2.2' are different")
      }
    }else{
      par2.2 <- par1.2
    }
  }else{
    if (is.null(par2.2)){
      stop("When 'var.equal' is set to FALSE, a value for 'par2.2' must be given")
    }
    if (par1.2==par2.2){
      stop("'var.equal' is set to FALSE, but par1.2==par2.2")
    }
  }
  costs.origin <- costs
  rho.origin <- rho
  par1.1.origin <- par1.1
  par2.1.origin <- par2.1
  L <- width/2
  if (par1.1 > par2.1){
    rho <- 1 - rho
    costs <- costs[, 2:1] # change c.t.pos <-> c.t.neg and c.f.pos <-> c.f.neg
    g <- par1.1; par1.1 <- par2.1; par2.1 <- g
    f <- par1.2; par1.2 <- par2.2; par2.2 <- f
  }
  # slope
  R <- slope(rho, costs)
  if(var.equal){
    num <- 1/2-par1.2^2*log(R)/(par2.1-par1.1)^2
    den <- 1/2+par1.2^2*log(R)/(par2.1-par1.1)^2
    epsilon <- abs(num/den)
    est.non.dis <- (qnorm(1-alpha/2)/L)^2*((log(R)/(par2.1-par1.1))^2*2*par1.2^4/(1+epsilon)
                                           +(1/2-par1.2^2*log(R)/(par2.1-par1.1)^2)^2*par1.2^2/epsilon
                                           +(1/2+par1.2^2*log(R)/(par2.1-par1.1)^2)^2*par1.2^2)
  }else{
    threshold <- expression((sigma2D*muND-sigma2ND*muD+sqrt(sigma2ND*sigma2D)*sqrt((muD-muND)^2+2*log(R*sqrt(sigma2D/sigma2ND))*(sigma2D-sigma2ND)))/(sigma2D-sigma2ND))
    where <- list(muND=par1.1, sigma2ND=par1.2^2, muD=par2.1, sigma2D=par2.2^2, R=R)
    a <- (as.numeric(attributes(eval(deriv(threshold, "muD"), where))$gradient))^2
    b <- (as.numeric(attributes(eval(deriv(threshold, "muND"), where))$gradient))^2
    c <- (as.numeric(attributes(eval(deriv(threshold, "sigma2D"), where))$gradient))^2
    d <- (as.numeric(attributes(eval(deriv(threshold, "sigma2ND"), where))$gradient))^2
    epsilon <- sqrt((a*par2.2^2+2*c*par2.2^4)/(b*par1.2^2+2*d*par1.2^4))
    est.non.dis <- (qnorm(1-alpha/2)/L)^2*(a*par2.2^2/epsilon+b*par1.2^2+2*c*par2.2^4/epsilon+2*d*par1.2^4)
  }
  est.dis <- epsilon*est.non.dis
  
  if (par1.1.origin > par2.1.origin){
    auxiliar <- est.non.dis
    est.non.dis <- est.dis
    est.dis <- auxiliar
    epsilon <- 1/epsilon
  }
  
  re <- list(ss2=est.dis, ss1=est.non.dis, epsilon=epsilon, width=width, alpha=alpha, costs=costs.origin, R=R, prev=rho.origin)
  class(re) <- "SS"
  return(re)
} 

#############################################################
## Print function for class "SS"
#############################################################
print.SS <- function(x, ...){
    cat("Optimum SS Ratio: ", x$epsilon)
    cat("\n\nSample size for")
    cat("\n  Diseased: ", x$ss2)
    cat("\n  Non-diseased: ", x$ss1)
    cat("\n")
    cat("\nParameters used")
    cat("\n  Significance Level: ", x$alpha)
    cat("\n  CI width: ", x$width)
    cat("\n  Disease prevalence:", x$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$costs)
    cat("\n  R:", x$R)
    cat("\n")
}


