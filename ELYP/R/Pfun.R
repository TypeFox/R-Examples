Pfun <- function(b1, b2, Mulam) {
####
#### The hazard ratio of (Z=1) over (Z=0) in the YP model is computed here.
#### b1 is the short term hazard ratio, b2 is the long term hazard ratio.
#### The Mulam is the value from myLLfun( ) and is defined by int g(t) dH(t)
#### with g(t) = I[t <= 2.0] etc.
####
alpha <- exp(-Mulam)
TrtCon <- 1/(alpha*exp(-b1) + (1-alpha)*exp(-b2))
return(TrtCon)
}
