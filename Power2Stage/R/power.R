#------------------------------------------------------------------------------
# functions for power calculation
# Author: dlabes
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# helper function to allow partial match of the method
.powerMethod <- function(method){
  meth <- tolower(method[1])
  if (method=="") method <- "exact"
  methods <- c("exact","owenq","noncentral","nct","shifted","central","mvt")
  #                         ^ = match at start
  meth <- methods[grep(paste("^",meth,sep=""), methods)]
  if (length(meth)==0) meth <- tolower(method)
  if (length(meth)>1)  meth <- "nct" # happens if only "n" is given
  
  return(meth)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks,
# does not vectorize propperly!
# to be used in sampleN.TOST avoiding overhead of redundant calculations
# in case of multiplicative model:
# diffm=log(null ratio), theta1=log(lower BE limit), theta2=log(upper BE limit)
# in case of additive model:
# diffm=1-null ratio, theta1=lower BE limit-1, theta2=upper BE limit -1
# Jan 2015: Interface changed to sem 
# so call it with sem= se*sqrt(bk/n) if balanced or se*sqrt(bkni*sum(1/n)) 
.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df)
{
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
  # and se=0!
  delta1 <- (diffm-ltheta1)/sem
  delta2 <- (diffm-ltheta2)/sem
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  # R is infinite in case of alpha=0.5 where tval is == 0
  R <- (delta1-delta2)*sqrt(df)/(2.*tval)
  # in case of se=0 it results: delta1=Inf, delta2=inf if diffm>ltheta2
  # Inf - Inf is NaN
  R[is.nan(R)] <- 0
  
  # if alpha>0.5 (very unusual!) t(1-alpha,df) is <0 and then R is negative 
  # i.e in OwensQ the upper integration limit is lower then the lower limit!
  # SAS OwenQ gives missings if b or a are negative!
  # On the other hand SAS Proc Power gives values which are seemingly calculated
  # with abs(R). 
  # Correct acc. to Fig. 1 given in K.Philips
  # "Power for Testing Multiple Instances of the Two One-Sided Tests Procedure"
  # The International Journal of Biostatistics: Vol. 5: Iss. 1, Article 15.
  # DOI: 10.2202/1557-4679.1169
  # should be R=Inf, i.e. unlimited integration with respect to sigma.
  # This gives the same values (within certain precision) as Ben's power.1TOST
  # aka power.TOST(..., method="mvt").
  # Can also be checked via function power.TOST.sim().

  R[R<=0] <- Inf
  # comment above out and write
  # R <- abs(R)
  # to check SAS Proc power values
  
  # to avoid numerical errors in OwensQ implementation
  if (min(df)>10000){
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    tval <- qnorm(1-alpha)
    p1   <- pnorm(tval-delta1)
    p2   <- pnorm(-tval-delta2)
    return(p2-p1)
  }
  if (min(df)>=5000 & min(df<=10000)) {
    # approximation via non-central t-distribution
    return(.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df))
  }
  
  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se OR n,df are vectors) 
  nel <- length(delta1)
  dl <- length(tval)
  p1 <- c(1:nel)	
  p2 <- p1
  for (i in seq_along(delta1)) {
    if (dl>1) {
      ddf <- df[i]; ttt <- tval[i]
    } else {
      ddf <- df[1]; ttt <- tval[1]
    }
    p1[i] <- OwensQ(ddf,  ttt, delta1[i], 0, R[i])
    p2[i] <- OwensQ(ddf, -ttt, delta2[i], 0, R[i])
  }
  #browser()
  pow <- p2-p1
  # due to numeric inaccuracies power < 0?
  pow[pow<0] <- 0
  return( pow )
}

#------------------------------------------------------------------------------
# Ben's implementation of power via integration of the bivariate t-distribution
# with correlation == 1, also exact
# does'nt vectorize in any respect!
.power.1TOST <- function(alpha, ltheta1, ltheta2, diffm, sem, df, setseed = TRUE)
{
  if (setseed) set.seed(123456)
  corr  <- matrix(1, ncol = 2, nrow = 2)
  
  tval  <- qt(1 - alpha, df)
  lower <- c(tval, -Inf)
  upper <- c(Inf, -tval)
  delta1 <- (diffm - ltheta1) / sem
  delta2 <- (diffm - ltheta2) / sem
  pow <- rep(0, times=length(delta1))
  # attempt to vectorize if ltheta0 OR se is a vector
  for(i in seq_along(delta1)){
    delta <- c(delta1[i], delta2[i])
    prob  <- pmvt(lower = lower, upper = upper, delta = delta, df = df, corr = corr,
                  algorithm = GenzBretz(maxpts=100000, abseps = 1e-05))#[1]
    # abseps=1e-6 gives often "Completion with error > abseps"
    # give a warning if attr(prob,"msg") not equal "Normal completion"?
    if(attr(prob, which="msg")!="Normal Completion")
      warning("pmvt returned message ", attr(prob, which="msg"), call.=FALSE)
    pow[i] <- prob[i]
  }
  pow
}

#------------------------------------------------------------------------------
# 'raw' approximate power function without any error checks, 
# approximation based on non-central t
# this vectorizes ok
.approx.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE, log.p = FALSE)
  
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
  # and sem=0!
  delta1 <- (diffm-ltheta1)/sem
  delta2 <- (diffm-ltheta2)/sem
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  
  pow <- pt(-tval, df, ncp=delta2)-pt(tval, df, ncp=delta1)
  pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
  
  return(pow)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks, 
# approximation based on central 'shifted' central t distribution
# according to Chow, Liu "Design and Analysis of Bioavailability ..."
# Chapter 9.6 and implemented in PASS 2008
# where does this all come from?
.approx2.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df)
{
	tval   <- qt(1 - alpha, df, lower.tail = TRUE)
	# 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
	# and se=0!
	delta1 <- (diffm-ltheta1)/sem
	delta2 <- (diffm-ltheta2)/sem
	# is this correct?
	delta1[is.nan(delta1)] <- 0
	delta2[is.nan(delta2)] <- 0
	
	pow <- pt(-tval-delta2, df) - pt(tval-delta1, df)
	pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
	
	return(pow)
}
#------------------------------------------------------------------------------
# function for merging the various power calculations
.calc.power <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df, method="exact")
{ 
  pow <- switch(
      method,
      exact=.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      owenq=.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      nct=  .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      noncentral=.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      shifted=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      central=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      mvt=.power.1TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      stop("Method '", method, "' unknown!\n", call.=TRUE)
  ) 
  return(pow)
}
