#-------------------------------------------------------------------------------
# author: dlabes
#------------------------------------------------------------------------------
#-- The functions of normal-, t-distributions and integrate() ------------------
#require(stats) #this is usually not necessary within a standard installation
# Owen's Q-function 
# a, b must be a scalar numeric
# nu, t and delta also, no vectors allowed
OwensQ <- function (nu, t, delta, a, b)
{
  if(missing(a)) a <- 0
  
	if (length(nu)>1 | length(t)>1 | length(delta)>1 | length(a)>1 | 
      length(b)>1) stop("Input must be scalars!")
  if (nu<1) stop("nu must be >=1!")
  
  if(a==b) return(0)
  
  # Observation: for 'really' large df (nu>2000) and large delta/b the  
	# density function is zero over nearly all its range! Q than returned 
	# sometimes falsly as =0! See documentation ?integrate for that.
	# example: OwensQ(3000,1.64,-10,0,300) = 1
	#          OwensQ(3000,1.64,-10,0,350) = 5.614476e-12 ! Correct is 1.
	# Idea: adapt upper and/or lower integration limit to account for that
	low <- a; up <- b
  
  # May 2011: shrink  interval depending on delta*b 
  #if (b>1E6) b <- Inf
	# in case of alpha=0.5 b is infinite
  # also in case of se=0 and diffm != ltheta1 or !=ltheta2
	if (is.finite(b)){
  	if (nu >= 1000 || abs(delta*b) > 30 || b>50){
      # all adaptions don't cover all extremal cases!
      
      # return OwensQOwens for cases with high b, delta*b
      # but only for df<400, else computation time may stuck
      # up to 400 a call is <= 2 msec
      if(nu<400 & a==0) return(OwensQOwen(nu,  t, delta, 0, b))
      
      # try to shorten the range via interval halving: Jul 2012, se below
      # but this doesn't work for Jiri's extremal cases
      # before we had checked the integrand at 500 points
      # this is from computational time better then interval halving
      # at least for high nu
      step <- (b-a)/499
      x <- a + (0:499) * step
      dens <- .Q.integrand(x, nu, t, delta)
      x <- x[dens!=0]
      if (length(x)>0) {
        low <- max(min(x) - step, a)
        up  <- min(max(x) + step, b)
      }
      # 
#       # upper integration limit via interval halving
#       ab <- b-a
#       x  <- b
#       dens <- .Q.integrand(x, nu, t, delta)
#       #cat("Upper search\n")
#       #cat("x=",x,"dens=",dens,"\n")
#       while (dens==0){
#         ab   <- ab*0.5
#         x    <- x - ab
#         dens <- .Q.integrand(x, nu, t, delta)
#         #cat("x=",x,"dens=",dens,"\n")
#         if (ab < 1E-10) break
#       }
#       if (dens>0) up <- min(x + ab, b) #else return(0)
#       #cat("up=",up,"\n")
#       # lower limit
#       ab   <- up - a
#       x    <- a
#       dens <- .Q.integrand(x, nu, t, delta)
#       #cat("Lower search\n")
#       #cat("x=",x,"dens=",dens,"\n")
#       while (dens==0){
#         ab   <- ab*0.5 
#         x    <- x + ab
#         dens <- .Q.integrand(x, nu, t, delta)
#         #cat("x=",x,"dens=",dens,"\n")
#         if (ab<1E-10) break
#       }
#       if (dens>0)  low <- max(x - ab, a) #else return(0)
#       #cat("low=",low,"\n")
    }  
  }
	# result of integrate() is a list, see ?integrate
	# .Machine$double.eps^.5 = 1.490116e-08 on my machine
	# MBESS uses .Machine$double.eps^0.25 = 0.0001220703 for both tolerances
	# seems it makes no difference
	Qintegral <- integrate(.Q.integrand, lower = low, upper = up, 
			         nu=nu, t=t, delta = delta, subdivisions = 10000, 
			         #rel.tol = .Machine$double.eps^0.5, 
               #abs.tol = .Machine$double.eps^0.5,
			         rel.tol = 1.e-9, abs.tol=1.e-9, stop.on.error = TRUE)[[1]]
	# error handling? How?
	if(a==0){
    # approx. via nct?
    # this seems not correct for all cases!
    # f.i. 
    #>OwensQ(nu=100,t=1.64, delta=1,a=0,b=1)
    # [1] 0.7361993
    # Warning message:
    #  In OwensQ(nu = 100, t = 1.64, delta = 1, a = 0, b = 1) :
    #  OwensQ = 3.597533e-81 due to numeric problems.
    #  Replaced with nct approx. = 0.7361993
    #> OwensQOwen(nu=100,t=1.64, delta=1,a=0,b=1)
    # [1] 8.326673e-17
    # but seems to work for sufficient high b
    # the question is what is sufficient
    #browser()
	  if (b>50){
      check <- pt(t, df=nu, ncp=delta)
      # debug print
      # cat("Check aginst nct approx.:", check, Qintegral,"\n")
      if(Qintegral<1e-9 & check>1e-7) {
        warning("OwensQ = ", format(Qintegral, digits=7),
                " due to numeric problems.",
                "\n  Replaced with nct approx. = ", format(check, digits=7))
        Qintegral <- check
      }  
    }
  }
	return(Qintegral)  
}
#-------------------------------------------------------------------------------
# Integrand of the definit integral in Owen's Q. Used in the call of integrate()
# Not useful alone, I think ? Leading . hides this function 
# function must give a vectorized answer in respect to x
.Q.integrand <- function(x, nu, t, delta)
{ #version without for - loop, it works without
	lnQconst <- -((nu/2.0)-1.0)*log(2.0) - lgamma(nu/2.)

# what if x<0? Should here not possible, but ...
# simple x^(nu-1) doesnt work for high nu because  = Inf 
# and then exp( -0.5*x^2 + lnQconst )*x^(nu-1) -> NaN
# (nu-1)*log(abs(x)) is NaN if nu=1, x=0! 0*(-Inf) -> NaN
  
  dens <- x  # assures that dens=0 if x=0
  x1   <- x[x!=0]
  dens[x!=0] <- sign(x1)^(nu-1) *
      pnorm( t*x1/sqrt(nu) - delta, mean = 0, sd = 1, log.p = FALSE) * 
      exp( (nu-1)*log(abs(x1)) - 0.5*x1^2 + lnQconst )
  # return    
	dens
}

# Test cases:
# Craig Zupke's observations:
# power.TOST(0.410,FALSE,-5.97,5.97,8.5448,1,14,"parallel",TRUE) #!old call
# power.TOST(0.410,FALSE,-5.97,5.97,8.5448,1,14,"parallel","exact")
# gave an error; high b/delta
# should give: 2.335633e-07

# Jul 2012: Helmuts observation
# n=4, CV=1E-5(=se) gives power=1 (delta1=24303.3, delta2=-38811.23, R=b=15283.88
#      CV=1E-6 gives power=0      (      243033          -388112.3   R  152838.8
#      CV=0    gives power=1             Inf              -Inf       Inf
# tval=2.919986
# for CV=1e-6: erroneous in versions pre 0.9-9. 2. call gave =0
# OwensQ(nu=2, t= 2.919986, delta= 243033,   0, 152838.8) ==0
# OwensQ(nu=2, t=-2.919986, delta=-388112.3, 0, 152838.8) ==1
# for CV=0 - no longer staistics
# OwensQ(nu=2, t=2.919986, delta=Inf, 0, Inf)  ==0
# OwensQ(nu=2, t=-2.919986, delta=-Inf, 0,Inf) ==1
# 
# Helmuts cases (ver 0.9-9) Jul 2012
# sampleN.TOST(theta0=1, CV=0.02, design="2x2", print=TRUE) # Ok
# next gave an error due to 0*-Inf in .Q.integrand()
# sampleN.TOST(theta0=1, CV=0.01, design="2x2", print=TRUE) 
#
# Jiri Hofmann's case(s) (V1.1-08, Jan 2014)
# power.TOST(CV=0.39, n=7000, theta0=1.24) gave power=0.3511889, correct
# power.TOST(CV=0.385, n=7000, theta0=1.24) gave power=0, correct is 0.3568784
# power.TOST(CV=0.38, n=7000, theta0=1.24) gave power=0, correct is 0.3627554
# power.TOST(CV=0.375, n=7000, theta0=1.24) gave power=0, correct is 0.3688277
# power.TOST(CV=0.37, n=7000, theta0=1.24) gave power=0.3751032, correct
#
# power.TOST(CV=0.37, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.375, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.38, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.385, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.39, n=7000, theta0=1.0) gave power=1, correct
