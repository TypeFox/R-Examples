# -----------------------------------------------------------------------------
# function to obtain exact critical values of bivariate nct
# Author B. Lang
# modified by D. Labes
# -----------------------------------------------------------------------------
cval_fixed <- function(se, df, ltheta1, ltheta2, alpha = 0.05) 
{
  tval <- qt(1 - alpha, df)
  crange <- c(-tval-0.1, tval+0.1)
  if (min(df) > 10000) return(tval)
  f <- function(c) {
    .power.TOST2(cval = c, alpha = alpha, ltheta1 = ltheta1, ltheta2 = ltheta2,
                 diffm = ltheta1, sem = se, df = df) - alpha
  }
  #browser()
  uniroot(f, interval = crange, extendInt = "yes")$root
}

# -----------------------------------------------------------------------------
# modified exact power function to use tval as critical values or 
# bivariate nct critical values
.power.TOST2 <- function(alpha=0.05, cval, ltheta1, ltheta2, diffm, sem, df)
{
  #browser()
  if(missing(cval)) tval <- qt(1 - alpha, df, lower.tail = TRUE) else tval <- cval
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2 and sem=0!
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
  # to check SAS Proc power values comment above out and write
  # R <- abs(R)
  
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
  pow <- p2-p1
  # due to numeric inaccuracies power < 0?
  # paranoia
  pow[pow<0] <- 0
  return( pow )
}
