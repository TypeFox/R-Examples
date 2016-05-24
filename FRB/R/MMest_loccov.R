MMest_loccov <- function(Y, control=MMcontrol(...), ...)
{
# computes multivariate location and shape (M)M-estimates with auxiliary S-scale 
# INPUT:
#   Y = data
# OUTPUT:
#   res$Mu : MM-location estimate
#   res$Gamma : MM-shape matrix
#   res$Sigma : MM-covariance
#   res$SMu : S-location estimate
#   res$SGamma : S-shape matrix
#   res$SSigma : S-shape matrix
#   res$scale : S-scale
#
# calls: Sest_loccov()
# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes scaled Tukey's biweight psi function with constant c, divided by x, for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------

convert <- function(control) {
# converts control list into an rrcov control object for MM-estimator

control.S=CovControlSest(eps = control$fastScontrols$convTol, 
maxiter = control$fastScontrols$maxIt,nsamp = control$fastScontrols$nsamp, 
trace = FALSE, method= "sfast")

control=CovControlMMest(bdp=control$bdp,eff=control$eff,sest=control.S,
tolSolve = control$convTol.MM, maxiter = control$maxIt.MM)
return(control)
}


# --------------------------------------------------------------------
# -                       main function                              -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
n <- nrow(Y)
m <- ncol(Y)


if (n<=m) stop("number of observations too small (should have n > m)")

MMests=CovMMest(Y,eff.shape=control$shapeEff,control=convert(control))

w <- scaledpsibiweight(sqrt(getDistance(MMests)),MMests@c1)


return(list( Mu = t(getCenter(MMests)), Gamma = getShape(MMests), 
scale = getDet(MMests)^(1/(2*m)), Sigma = getCov(MMests), 
c1=MMests@c1, 
SMu=t(getCenter(MMests@sest)),SGamma = getShape(MMests@sest),
SSigma = getCov(MMests@sest),c0=MMests@sest@cc,
b=MMests@sest@kp, w=w, outFlag=(!getFlag(MMests))))
}


