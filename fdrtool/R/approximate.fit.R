### approximate.fit.R  (2007-10-19)
###
###     First Guess of Null Model Parameters
###
### Copyright 2007 Korbinian Strimmer 
###
###
### This file is part of the `fdrtool' library for R and related languages.
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


# first guess of null model parameters
approximate.fit = function(x, statistic=c("normal", "correlation", "pvalue", "studentt"), lambda=seq(0,0.9,0.05))
{
  statistic <- match.arg(statistic)
  nm = get.nullmodel(statistic)

  param = iqr.fit(x, statistic)  # approximate estimate of scale parameter
  pval = nm$get.pval(x, param)   # compute corresponding p-values

  # estimate eta0
  eta0 = pval.estimate.eta0(pval, method="quantile", q=0.1, diagnostic.plot=FALSE, lambda=lambda) 

  return(list(eta0=eta0, param=param))
}


# find robust estimate of scale by fitting IQR
iqr.fit = function(x, statistic=c("normal", "correlation", "pvalue", "studentt"))
{
  statistic <- match.arg(statistic)
  nm = get.nullmodel(statistic)
  
  # observed interquantile range
  iqr.obs = as.double(diff(quantile(x, probs=c(.25, .75))))
  
  if (statistic=="pvalue")
  {
    param = NULL
  }
  else if (statistic=="normal")
  {
    param = iqr.obs/(2*qnorm(.75))
  }
  else
  {
    mfun = function(param) return( (nm$iqr(param)-iqr.obs)^2 )
    supp = nm$get.support()
    param = optimize(mfun, supp)$minimum  
  }

  return(param)
}



