
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
# FUNCTION:                  EMPIRICAL COPULAE PROBABILIY:
#  pempiricalCopula           Computes empirical copula probability
# FUNCTION:                  EMPIRICAL COPULAE DENSITY:
#  dempiricalCopula           Computes empirical copula density
# FUNCTION:                  DEBYE FUNCTION:
#  .Debye                     Returns the value of the Debye function of order k
#  .Debye1                    Returns the value of the Debye function of order 1
# FUNCTION:
#  .pmoCopula
#  .dmoCopula
################################################################################


################################################################################
# FUNCTION:                  EMPIRICAL COPULAE PROBABILIY:
#  pempiricalCopula           Computes empirical copula probability


pempiricalCopula <-
    function(u, v, N = 10)
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Computes the empirical copula probability
    
    # Source:
    #   bouye02a.pdf
    
    # FUNCTION:
    
    # Settings:
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 1]
        u = u[, 2]
    }
    
    # Probability:
    p = q = (0:N)/N
    h = matrix(rep(0, times = (N+1)^2), N+1)
    for ( i in (0:N) ) {
        for ( j in (0:N) ) {
            z = Heaviside(u, p[i+1]) + Heaviside(v, q[j+1])
            h[j+1, i+1] = length(z[z == 0])
        }
    }
    h = h/length(u)
    
    # Return Value:
    list(x = p, y = q, z = h)
}


################################################################################
# FUNCTION:                  EMPIRICAL COPULAE DENSITY:
#  dempiricalCopula           Computes empirical copula density


dempiricalCopula <-
    function(u, v, N = 10)
{   
    # A function implemented by Diethelm Wuertz

    # Description
    #   Computes the empirical copula probability
    
    # Source:
    #   bouye02a.pdf
    
    # FUNCTION:
    
    # Settings:
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 1]
        u = u[, 2]
    }
    
    # Probability:
    ans = pempiricalCopula(u, v, N)
    X = ans$x
    Y = ans$y
    C = ans$z
    
    # Density:
    M = N+1
    x = X[-1] - diff(X)/2
    y = Y[-1] - diff(Y)/2
    c = C[-1,-1]+C[-M,-M]-C[-1,-M]-C[-M,-1]
    
    # Return Value:
    list(x = x, y = y, z = c)
}


################################################################################
# FUNCTION:                  DEBYE FUNCTION:
#  .Debye                     Returns the value of the Debye function of order k
#  .Debye1


.Debye <- 
    function(x, k = 1) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the value of the Debye function of order k
    
    # Arguments:
    #   x - a numeric value or vector
    #   k - the order of the Debye function, a positive integer value   
    
    # FUNCTION:
    
    # Check:
    if (!is.integer(k) | k <= 0) 
        stop("k must be a positive integer")
    
    # Loop:
    D = NULL
    error = NULL
    for ( i in 1:length(x) ) {
        nextD = .Debye1(x[i],k) 
        D = c(D, nextD[[1]])
        error = c(error, nextD[[2]])
    }
    
    # Add error attribute:
    attr(D, "error") = error
    
    # Return Value:
    D
}


# ------------------------------------------------------------------------------


.Debye1 <-
    function(x, k = 1)
{   
    # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Function to be integrated:
    d = function(x, lambda) { x^lambda / ( exp(x) - 1 ) }   
    
    # Integrate:
    u = abs(x)
    if (x == 0) {
        D = 1
        error = 0   
    } else {
        ans = integrate(f = d, lower = 0, upper  = u, lambda = k)
        D = k * ans[[1]] / u^k
        error = ans[[2]]
    }
    if (x < 0) {
        D = D + k*u/(k+1)
    }
    
    # Return Value:
    list(D = D, error = error)
}


################################################################################
# FUNCTION:
#  .pmoCopula
#  .dmoCopula


.pmoCopula <- 
    function(u = 0.5, v = u, alpha = NULL)
{
    if (is.null(alpha)) alpha = c(0.5, 0.5)
    alpha1 = alpha[1]
    alpha2 = alpha[2]
    
    U = u^(1-alpha1) * v
    V = u * v^(1-alpha2)
    UV = cbind(U,V)
    apply(UV, 1, max)
}   


# ------------------------------------------------------------------------------


.dmoCopula <- 
    function(u = 0.5, v = u, alpha = NULL)
{
    if (is.null(alpha)) alpha = c(0.5, 0.5)
    alpha1 = alpha[1]
    alpha2 = alpha[2]
    
    U = u^(1-alpha1) * v
    V = u * v^(1-alpha2)
    UV = cbind(U,V)
    apply(UV, 1, max)
}   


################################################################################

