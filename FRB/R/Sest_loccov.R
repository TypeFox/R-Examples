Sest_loccov <- function(Y, bdp=.5, control=Scontrol(...),...)
{
# Computes S-estimates by using the 
# fast S algorithm for multivariate location/covariance estimation
# implemented in the rrcov package 
# INPUT:
#   Y : response matrix (n x m)
#   bdp : breakdown point (<= 0.5)
# OUTPUT:
#   res$Mu : estimate of regression coefficients (or location vector)
#   res$Gamma : estimate of shape matrix 
#   res$scale : estimate of scale

#---------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c, divided by x, for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------

vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------

reconvec <- function(vec,ncol) {
# reconstructs vecop'd matrix

lcol <- length(vec)/ncol
rec <- matrix(0,lcol,ncol)
for (i in 1:ncol)
    rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]

return(rec)
}

# --------------------------------------------------------------------

convert <- function(control) {
# converts control list into an rrcov control object for S-estimator

control=CovControlSest(eps = control$convTol, maxiter = control$maxIt,
        nsamp = control$nsamp, trace = FALSE, method= "sfast")
return(control)
}


# --------------------------------------------------------------------
# -                       main function                              -
# --------------------------------------------------------------------

#set.seed(10) # seed can be set, but be careful in simulations then...

Y <- as.matrix(Y)
n <- nrow(Y)
m <- ncol(Y)

if (n<=m) stop("number of observations too small (should have n > m)")

Sests <- CovSest(Y,bdp=bdp,control=convert(control))
#Smu <- getCenter(Sests)
#SSigma <- getCov(Sests)
w <- scaledpsibiweight(sqrt(getDistance(Sests)),Sests@cc)
return(list( Mu = t(getCenter(Sests)), Gamma = getShape(Sests), 
scale = getDet(Sests)^(1/(2*m)), Sigma = getCov(Sests), 
c=Sests@cc, b=Sests@kp, w=w, outFlag=(!getFlag(Sests))))
}


