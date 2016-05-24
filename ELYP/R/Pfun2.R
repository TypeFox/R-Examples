Pfun2 <- function(b1, b2, a, X, Mulam) {
####
#### The hazard ratio of (Z=1) over (Z=0) in the YP model is computed here.
#### b1 is the short term hazard ratio, b2 is the long term hazard ratio.
#### The Mulam is the value from myLLfun( ) and is defined by int g(t) dH(t)
#### with g(t) = I[t <= 2.0] etc.
#### aX is the part that is proportional hazard. a is parameter X is covariate.
#### This is try to reduce the dim. Since aX is always 1 dim.  aX = a^t X
####   this function is used by findU4( )
####
alpha <- exp(-Mulam)
TrtCon <- exp( a*X )/(alpha*exp(-b1) + (1-alpha)*exp(-b2))
return(TrtCon)
}
