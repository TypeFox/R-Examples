
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                    DESCRIPTION:
#  ghMean                       Returns true GH mean
#  ghVar                        Returns true GH variance
#  ghSkew                       Returns true GH skewness
#  ghKurt                       Returns true GH kurtosis
# FUNCTION:                    DESCRIPTION:
#  ghMoments                    Returns true GH moments
# FUNCTION:                    UTILITY FUNCTION:
#  .aRecursionGH                Computes the moment coefficients a recursively
#  .besselZ                     Computes Bessel/Power Function ratio
# FUNCTION:                    MOMENTS ABOUT MU:
#  .ghMuMoments                 Computes mu-moments from formula
#  .ghMuMomentsIntegrated       Computes mu-moments by integration
#  .checkGHMuMoments            Checks mu-moments
# FUNCTION:                    MOMENTS ABOUT ZERO:
#  .ghRawMoments                Computes raw moments from formula
#  .ghRawMomentsIntegrated      Computes raw moments by Integration
#  .checkGHRawMoments           Checks raw moments
# FUNCTION:                    MOMENTS ABOUT MEAN:
#  .ghCentralMoments            Computes central moments from formula
################################################################################


ghMean <-
function(alpha=1, beta=0, delta=1, mu=0, lambda=-1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    zeta = delta * sqrt(alpha*alpha-beta*beta)
    mean = mu + beta * delta^2 * .kappaGH(zeta, lambda)
    mean
}


# ------------------------------------------------------------------------------


ghVar <- 
function(alpha=1, beta=0, delta=1, mu=0, lambda=-1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    var = .ghCentralMoments(k=2, alpha, beta, delta, mu, lambda)[[1]]
    var
}


# ------------------------------------------------------------------------------


ghSkew <- 
function(alpha=1, beta=0, delta=1, mu=0, lambda=-1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Moments:
    k2 = .ghCentralMoments(k=2, alpha, beta, delta, mu, lambda)[[1]] 
    k3 = .ghCentralMoments(k=3, alpha, beta, delta, mu, lambda)[[1]]
    
    # Return Value:
    skew = k3/(k2^(3/2))     
    skew          
}


# ------------------------------------------------------------------------------


ghKurt <- 
function(alpha=1, beta=0, delta=1, mu=0, lambda=-1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Moments:
    k2 = .ghCentralMoments(k=2, alpha, beta, delta, mu, lambda)[[1]]
    k4 = .ghCentralMoments(k=4, alpha, beta, delta, mu, lambda)[[1]]

    # Return Value:
    kurt = k4/k2^2 - 3 
    kurt
}


################################################################################


ghMoments <-
function(order, type = c("raw", "central", "mu"),
    alpha=1, beta=0, delta=1, mu=0, lambda=-1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns moments of the GH distribution
    
    # FUNCTION:
    
    # Settings:
    type = match.arg(type)
    
    # Moments:
    if (type == "raw") {
        ans = .ghRawMoments(order, alpha, beta, delta, mu, lambda)
    } else if (type == "central") {
        ans = .ghCentralMoments(order, alpha, beta, delta, mu, lambda)
    } else if (type == "mu") {
        ans = .ghMuMoments(order, alpha, beta, delta, mu, lambda)  
    }
    
    # Return Value:
    ans   
}


################################################################################


.aRecursionGH <- 
function(k = 12, trace = FALSE)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the moment coefficients recursively
    
    # Example:
    #   .aRecursionGH()
    #     [,1] [,2] [,3] [,4] [,5]  [,6]   [,7]   [,8]   [,9]  [,10] [,11] [,12] 
    # [1,]  1    0    0    0    0     0      0      0      0      0     0     0 
    # [2,]  1    1    0    0    0     0      0      0      0      0     0     0 
    # [3,]  0    3    1    0    0     0      0      0      0      0     0     0 
    # [4,]  0    3    6    1    0     0      0      0      0      0     0     0 
    # [5,]  0    0   15   10    1     0      0      0      0      0     0     0 
    # [6,]  0    0   15   45   15     1      0      0      0      0     0     0 
    # [7,]  0    0    0  105  105    21      1      0      0      0     0     0 
    # [8,]  0    0    0  105  420   210     28      1      0      0     0     0 
    # [9,]  0    0    0    0  945  1260    378     36      1      0     0     0 
    #[10,]  0    0    0    0  945  4725   3150    630     45      1     0     0 
    #[11,]  0    0    0    0    0 10395  17325   6930    990     55     1     0 
    #[12,]  0    0    0    0    0 10395  62370  51975  13860   1485    66     1 
    
    # FUNCTION:
    
    # Setting Start Values:
    a = matrix(rep(0, times = k*k), ncol = k) 
    a[1, 1] = 1
    if (k > 1) a[2, 1] = 1
    
    # Compute all Cofficients by Recursion:
    if (k > 1) {
        for (d in 2:k) {
            for (l in 2:d) {
                a[d,l] = a[d-1, l-1] + a[d-1, l]*(2*l+1-d)
            }
        }
    }
    rownames(a) = paste("k=", 1:k, sep = "")
    colnames(a) = paste("l=", 1:k, sep = "")
    
    # Trace:
    if (trace) {
        cat("\n")
        print(a)
        cat("\n")
    }
    for (K in 1:k) {
        L = trunc((K+1)/2):K  
        M = 2*L-K
        s1 = as.character(a[K, L])
        s2 = paste("delta^", 2*L, sep = "")
        s3 = paste("beta^", M, sep = "")
        s4 = paste("Z_lambda+", L, sep = "" )
        if (trace) {
            cat("k =", K, "\n")
            print(s1)
            print(s2)
            print(s3)
            print(s4)
            cat("\n")
        }
    }
    if (trace) cat("\n")
    
    # Return Value:
    list(a = a[k, L], delta = 2*L, beta = M, besselOffset = L) 
}
    
    
# ------------------------------------------------------------------------------
 

.besselZ <- 
function(x, lambda) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes Bessel/Power Function ratio
    
    # Ratio:
    ratio = besselK(x, lambda)/x^lambda
    
    # Return Value:
    ratio
}


################################################################################

    
.ghMuMoments <-
function(k = 4, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes mu moments from formula
    
    # Note:
    #   mu is not used.
    
    # FUNCTION:
    
    # Check parameters:
    if (lambda >= 0) stopifnot(abs(beta) < alpha)
    if (lambda <  0) stopifnot(abs(beta) <= alpha)
    if (lambda >  0) stopifnot(delta >= 0)
    if (lambda <= 0) stopifnot(delta >  0)
    
    zeta = delta*sqrt(alpha^2-beta^2)
    cLambda = 1/.besselZ(zeta, lambda)
    
    a = .aRecursionGH(k)
    coeff = a$a
    deltaPow = a$delta
    betaPow = a$beta
    besselIndex = lambda + a$besselOffset
        
    M = coeff * delta^deltaPow * beta^betaPow * .besselZ(zeta, besselIndex)
    ans = cLambda * sum(M)
    attr(ans, "M") <- M
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.ghMuMomentsIntegrated <-
function(k = 4, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes mu moments by integration
    
    # FUNCTION:
    
    # Check parameters:
    if (lambda >= 0) stopifnot(abs(beta) < alpha)
    if (lambda <  0) stopifnot(abs(beta) <= alpha)
    if (lambda >  0) stopifnot(delta >= 0)
    if (lambda <= 0) stopifnot(delta >  0)
      
    mgh <- function(x, k, alpha, beta, delta, lambda) {        
        x^k * dgh(x, alpha, beta, delta, mu=0, lambda) }
    
    M = integrate(mgh, -Inf, Inf, k = k, alpha = alpha,
        beta = beta, delta = delta, lambda = lambda,
        subdivisions = 1000, rel.tol = .Machine$double.eps^0.5)[[1]]
    
    # Return Value:
    M
}


# ------------------------------------------------------------------------------


.checkGHMuMoments <-
function(K = 10)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Checks mu moments
    
    # Example:
    #   .checkGHMuMoments()
    
    # FUNCTION:
    
    # Compute Raw Moments:
    mu = NULL
    for (k in 1:K) {
        mu = c(mu, 
            .ghMuMoments(k, alpha = 1.1, beta = 0.1, delta = 0.9, mu = 0.1))
    }
    names(mu) = 1:K
    int = NULL
    for (k in 1:K) {
        int = c(int, .ghMuMomentsIntegrated(k, 
            alpha = 1.1, beta = 0.1, delta = 0.9, mu = 0.1))
    }
    names(int) = NULL
    
    # Return Value:
    rbind(mu, int)
}


# ------------------------------------------------------------------------------


.ghRawMoments <-
function(k = 4, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes mu-moments from formula
    
    # FUNCTION:
    
    # Settings:
    binomCoeff = choose(k, 0:k)
    
    # Compute Mu Moments:
    M = 1
    for (l in 1:k) M = c(M, .ghMuMoments(l, alpha, beta, delta, mu, lambda))
    muPower = mu^(k:0)
    muM = binomCoeff * muPower * M
    
    # Return Value:
    sum(muM)  
}


# ------------------------------------------------------------------------------


.ghRawMomentsIntegrated <- 
function(k = 4, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes mu-moments by integration

    # FUNCTION:
    
    # Check parameters:
    if (lambda >= 0) stopifnot(abs(beta) < alpha)
    if (lambda <  0) stopifnot(abs(beta) <= alpha)
    if (lambda >  0) stopifnot(delta >= 0)
    if (lambda <= 0) stopifnot(delta >  0)
        
    # Compute Mu Moments by Integration:
    mgh <- function(x, k, alpha, beta, delta, mu, lambda) {
        x^k * dgh(x, alpha, beta, delta, mu, lambda)  }
    M = integrate(mgh, -Inf, Inf, k = k, alpha = alpha,
        beta = beta, delta = delta, mu = mu, lambda = lambda,
        subdivisions = 1000, rel.tol = .Machine$double.eps^0.5)[[1]]
    
    # Return Value:
    M
}


# ------------------------------------------------------------------------------


.checkGHRawMoments <-
function(K = 10)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Checks raw-moments 
    
    # Example:
    #   .checkGHRawMoments()
    
    # FUNCTION:
    
    raw = NULL
    for (k in 1:K) {
        raw = c(raw, .ghRawMoments(k, 
            alpha = 1.1, beta = 0.1, delta = 0.9, mu = 0.1))
    }
    names(raw) = 1:K
    int = NULL
    for (k in 1:K) {
        int = c(int, .ghRawMomentsIntegrated(k, 
            alpha = 1.1, beta = 0.1, delta = 0.9, mu = 0.1))
    }
    names(int) = NULL
    
    # Return Value:
    rbind(raw,int)
}


################################################################################
    
  
.ghCentralMoments <-
function (k = 4, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2) 
{
    zeta = delta * sqrt(alpha*alpha-beta*beta)
    mean = mu + beta * delta^2 * .kappaGH(zeta, lambda)
    M = 1
    for (i in 1:k)
        M = c(M, .ghMuMoments(i, alpha, beta, delta, mu, lambda)) 

    binomCoeff <- choose(k, 0:k)
    centralPower <- (mu - mean)^(k:0)
    centralM <- binomCoeff * centralPower * M 
    sum(centralM)
}


################################################################################

