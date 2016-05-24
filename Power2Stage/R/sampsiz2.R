# -------------------------------------------------------------------------
# sample size function without all the overhead
# for 2x2 crossover or 2-group parallel
# power via nct or shifted t approximation or exact
#
# Version with vectorize sample size w.r.t mse and/or ltheta0
#
# author D. Labes Jan 2015
# -------------------------------------------------------------------------
# library(PowerTOST)
# source("./R/sampsiz_n0.R")
# is also used for parallel groups with bk=4

.sampleN2 <- function(alpha=0.05, targetpower=0.8, ltheta0, ltheta1=log(0.8), 
                      ltheta2=log(1.25), mse, method="nct", bk=2)
{

  # se and ltheta0/diffm must have the same length to vectorize propperly!
  if (length(mse)==1)     mse <- rep(mse, times=length(ltheta0))
  if (length(ltheta0)==1) ltheta0 <- rep(ltheta0, times=length(mse))
  
  # return 'Inf' if ltheta0 not between or very near to ltheta1, ltheta2
  ns <- ifelse((ltheta0-ltheta1)<1.25e-5 | (ltheta2-ltheta0)<1.25e-5, Inf, 0)

  # design characteristics for 2-group parallel and 2x2 crossover design
  steps <- 2     # stepsize for sample size search
  nmin  <- 4     # minimum n
  
  se    <- sqrt(mse[is.finite(ns)])
  diffm <- ltheta0[is.finite(ns)]
  
  # start value from large sample approx. (hidden func.)
  # Jan 2015 changed to modified Zhang's formula
  # gives at least for 2x2 the best estimate (max diff to n: +-4)
  n <- .sampleN0_3(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps, bk)
  n <- ifelse(n<nmin, nmin, n)
  
  # method=="ls" is not used yet, in power.2stage.ssr() the 'original' ls approx
  # is used exactly as given in the paper of Golkowski et al.
  if(method=="ls") return(n)
  
  # degrees of freedom as expression
  # n-2 for 2x2 crossover and 2-group parallel design
  dfe <- parse(text="n-2", srcfile=NULL)
  # or should that read n-3? see Kieser/Rauch
  #dfe <- parse(text="n-3", srcfile=NULL)
  
  df   <- eval(dfe)
  pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df, method)
  iter <- rep(0, times=length(se)) 
  imax <- 50
  # imax =max. iter >50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max. of 2-3 steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
#  down <- FALSE; up <- FALSE
  index <- pow>targetpower & n>nmin
  while (any(index)) {
#    down <- TRUE
    n[index]    <- n[index]-steps     # step down if start power is to high
    iter[index] <- iter[index]+1
    df <- eval(dfe)
    pow[index]  <- .calc.power(alpha, ltheta1, ltheta2, diffm[index], 
                               sem=se[index]*sqrt(bk/n[index]), df[index],
                               method)
    index <- pow>targetpower & n>nmin & iter<=imax
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }

  # --- loop until power >= target power
  index <- pow<targetpower
  while (any(index)) {
#    up   <- TRUE; down <- FALSE
    n[index] <- n[index]+steps
    iter[index] <- iter[index]+1
    df <- eval(dfe)
    pow[index]  <- .calc.power(alpha, ltheta1, ltheta2, diffm[index], 
                               sem=se[index]*sqrt(bk/n[index]), df[index], 
                               method)
    index <- pow<targetpower & iter<=imax
  }
#   nlast <- n
#   if ((up & pow<targetpower) | (down & pow>targetpower) ) {
#     n <- NA
#   }
  # combine the Inf and n
  ns[ns==0] <- n
  n <- ns
  #  return only n here
  return(n)
  
} # end of function
