
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA






################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterArma             Filter function for Arma, Ar or Ma process
#                          
################################################################################


.filterArma <- function(
  data, size.data,
  m,n, 
  mu, a, b)
{  
    # Arguments
    # data: vector with data
    # m,n: model order arma(m,n)
    # mu,a,b: model parameters. 'a' is the autorregressive coefficient 
    # and 'b' is the moving average coeficient.
    # REMARKS: This function filters, pure 'arma' models, but 
    # it can deals with other models, such as 'ar' and 'ma' models. 
    # All the model parameters must be specified, regardless if they
    # exist. For example:
    # for ar(m) m >= 1 make n = 1 and b = 0.
    # for ma(n) n >= 1 make m = 1 and a = 0.
    
    # Return
    # the residuals 'z'
    
    # FUNCTION:
    
    # error treatment of input parameters
    if( n < 1 || n < 1 || length(a) != m || length(b) != n)
        stop("One or more of these conditions were true:
           n < 1 || n < 1 || length(a) != m || length(b) != n")
    
    # Initial declaration of variables
    x.init <- rep(0,m)
    z.init <- rep(0,n) 
    
    N <- size.data
    x.ZeroMean <- data - mu
    x2 = 0
    V <- c(x.init,x.ZeroMean[1:(N-1)])
    for( i in 1:m)
    {
      x1 <- -a[i]*V
      x2 = x2 +  x1[(m-(i-1)):(m+N-i)]
    }
    x2 <- x.ZeroMean + x2
    
    z <- filter(x2, filter = -b,
                method = "recursive", init = z.init)        
    
    if(length(z) != N)
        stop("Error in filtering function. length(z) != N")    
    
    # return
    return(z)
}


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterAparch          Filtering Aparch process and its particular cases
#                          
################################################################################

.filterAparch <- function(
  data,init.Value = NULL,
  p,q, 
  mu, omega, alpha, beta, gamma, delta)
{
  
    # Arguments
    # data: vector with data
    # p,q: model order. aparch(p,q)
    # mu,omega,alpha,beta,gamma,delta: model parameters
    # REMARKS: This function filters, pure 'aparch' models, but 
    # it can deals with other models, such as garch, arch. 
    # All the model parameters must be specified, regardless if they
    # exist. For example:
    # for garch(p,q) with p,q >= 1 make gamma = 0.
    # for arch(p) p >= 1 make q = 1, beta = 0 and gamma = 0.
    
    # Return
    # the series 'z' and 'h'
  
  
    # input 
    # data: vector with data
  
    # error treatment of input parameters
    if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
        || length(beta) != q)
        stop("One or more of these conditions were true:
          p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
          || length(beta) != q")
    
    # Initial declaration of variables
    pq = max(p,q)
    z = (data-mu)
    N = length(data)
    Mean.z = mean(abs(z)^delta)
    
    # initializing the time series
    if(is.null(init.Value)) {
        edelta.init <- rep(Mean.z,p)
        h.init <- rep(Mean.z,q)
    } else {
        edelta.init <- rep(init.Value,p)
        h.init <- rep(init.Value,q)  
    }  
    
    edeltat = 0
    for( i in 1:p)
    {
      edelta <- alpha[i]*(c(edelta.init,((abs(z)-gamma[i]*z)^delta)[1:(N-1)]))
      edeltat = edeltat +  edelta[(p-(i-1)):(p+N-i)]
    }
    edeltat = omega + edeltat
    
    h <- filter(edeltat, filter = beta,
                method = "recursive", init = h.init)
    
    if(length(z) != length(h))
      stop("Error in filtering function. length(z) != length(h)")    
    hh <- (abs(h))^(1/delta)
    
    # return
    cbind(z,hh)
}

################################################################################









################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterGarch11Fit      Filter function for Garch(1,1) process. This code 
#                          snipet was taken from the Wurtz et al. (2006) paper
#                          
################################################################################

.filterGarch11Fit <- function(data, parm)
{
  mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
  z = (data-mu); Mean = mean(z^2)
  
  # Use Filter Representation:
  e = omega + alpha * c(Mean, z[-length(data)]^2)
  h = filter(e, beta, "r", init = Mean)
  print(h[1])
  hh = sqrt(abs(h))
  cbind(z,hh)  
}

################################################################################













################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterAparchForLoop   Conditional Variance filtering - 
#                          for loop - Wuertz et al. (2006) 
#                          
################################################################################

.filterAparchForLoop <- function(
    data,h.init = 0.1,
    p,q, 
    mu, omega, alpha, beta, gamma, delta)
{
  
    # Return
    # the series 'z' and 'hh'
    
    
    # input 
    # data: vector with data
    
    # error treatment of input parameters
    if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
        || length(beta) != q)
      stop("One or more of these conditions were true:
            p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
            || length(beta) != q")
  
    # Initial declaration of variables
    pq = max(p,q)
    
    z = (data-mu)
    N = length(data)
    Mean.z = mean(abs(z)^delta)
    h = rep(h.init, pq)
    
    # Calculate h[(pq+1):N] recursively
    for (i in (pq+1):N )
    {
        ed = 0
        for (j in 1:p)
        {        
            ed = ed+alpha[j]*(abs(z[i-j])-gamma[j]*z[i-j])^delta
        }
        h[i] = omega + ed + sum(beta*h[i-(1:q)])
    }
    if(length(z) != length(h))
      stop("Error in filtering function. length(z) != length(h)")  
    
    hh <- (abs(h))^(1/delta)
    
    # return
    cbind(z,hh)    
}

################################################################################




