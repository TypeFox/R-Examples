### entropy.shrink.R  (2013-07-16)
###
###    Shrinkage estimators of entropy, mutual information
###    and related quantities
###
### Copyright 2008-13 Korbinian Strimmer
###
###
### This file is part of the `entropy' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


# shrinkage estimate of entropy 

# y:  a vector of counts (may include zeros)
entropy.shrink = function(y, lambda.freqs, unit=c("log", "log2", "log10"), verbose=TRUE)
{
  f = freqs.shrink(y, lambda.freqs=lambda.freqs, verbose=verbose)
  h = entropy.plugin(f, unit=unit)
  attr(h, "lambda.freqs") = attr(f, "lambda.freqs") # shrinkage intensity

  return( h )
}

freqs.shrink = function (y, lambda.freqs, verbose = TRUE) 
{
    target = 1/length(y) # uniform target (note length() works also for matrices)
    n = sum(y)
    u = y/n

    if (missing(lambda.freqs))
    {
      if (n==1 || n==0)
      {
        lambda.freqs = 1
      }
      else
      {
        lambda.freqs = get.lambda.shrink(n, u, target, verbose)
      }

    }
    else
    {
      if (verbose)
      {
        cat(paste("Specified shrinkage intensity lambda.freq (frequencies):", 
           round(lambda.freqs, 4)) , "\n")
      }
     
    }
    u.shrink = lambda.freqs * target + (1 - lambda.freqs) * u

    attr(u.shrink, "lambda.freqs") = lambda.freqs
    
    return(u.shrink)
}


# shrinkage estimation of mutual information
mi.shrink = function(y2d, lambda.freqs, unit=c("log", "log2", "log10"), verbose=TRUE)
{
  f2d = freqs.shrink(y2d, lambda.freqs=lambda.freqs, verbose=verbose)
  mi = mi.plugin(f2d, unit=unit)
  attr(mi, "lambda.freqs") = attr(f2d, "lambda.freqs") # shrinkage intensity

  return( mi )
}

# shrinkage estimation of chi-squared of independence
chi2indep.shrink = function(y2d, lambda.freqs, unit=c("log", "log2", "log10"), verbose=TRUE)
{
  f2d = freqs.shrink(y2d, lambda.freqs=lambda.freqs, verbose=verbose)
  chi2 = chi2indep.plugin(f2d, unit=unit)
  attr(chi2, "lambda.freqs") = attr(f2d, "lambda.freqs") # shrinkage intensity

  return( chi2 )
}


# shrinkage estimation of chi-squared statistic
chi2.shrink = function(y1, y2, lambda.freqs1, lambda.freqs2,
                       unit=c("log", "log2", "log10"), verbose=TRUE)
{
  f1 = freqs.shrink(y1, lambda.freqs=lambda.freqs1, verbose=verbose)
  f2 = freqs.shrink(y2, lambda.freqs=lambda.freqs2, verbose=verbose)
  chi2 = chi2.plugin(f1, f2, unit=unit)
  attr(chi2, "lambda.freqs1") = attr(f1, "lambda.freqs") # shrinkage intensity 1
  attr(chi2, "lambda.freqs2") = attr(f2, "lambda.freqs") # shrinkage intensity 2

  return( chi2 )
}

# shrinkage estimation of KL divergence
KL.shrink = function(y1, y2, lambda.freqs1, lambda.freqs2,
                     unit=c("log", "log2", "log10"), verbose=TRUE)
{
  f1 = freqs.shrink(y1, lambda.freqs=lambda.freqs1, verbose=verbose)
  f2 = freqs.shrink(y2, lambda.freqs=lambda.freqs2, verbose=verbose)
  KL = KL.plugin(f1, f2, unit=unit)
  attr(KL, "lambda.freqs1") = attr(f1, "lambda.freqs") # shrinkage intensity 1
  attr(KL, "lambda.freqs2") = attr(f2, "lambda.freqs") # shrinkage intensity 2

  return( KL )
}




## private function

get.lambda.shrink = function(n, u, t, verbose)
{
  # *unbiased* estimator of variance of u
  varu = u*(1-u)/(n-1)
  
  # misspecification
  msp = sum( (u-t)^2 )

  # estimate shrinkage intensity  
  if (msp == 0)
  {
    #warning("Overshrinkage")
    lambda = 1
  }
  else
    lambda = sum( varu ) / msp
  
  if (lambda > 1)
  {
    lambda = 1 # truncate at 1
    #warning("Overshrinkage")
  }
  
  if (lambda < 0)
  {
    lambda = 0
    #warning("Undershrinkage")
  }
  
  if (verbose)
  {
    cat(paste("Estimating optimal shrinkage intensity lambda.freq (frequencies):", 
      round(lambda, 4)) , "\n")
  }

  return(lambda)
}




