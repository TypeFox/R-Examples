### hc.score.R  (2012-07-25)
###
###     Compute empirical higher criticism score and threshold from p-values
###
### Copyright 2010-2012 Bernd Klaus and Korbinian Strimmer 
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



# compute empirical HC scores from p-values
hc.score = function(pval)
{
  if (max(pval) > 1 | min(pval) < 0) 
      stop("input p-values must all be in the range 0 to 1!")

  d = length(pval)
  F = rank(pval, ties.method="max")/d
  
  v = F*(1-F)/d # variance
  v[v==0] = min( v[v > 0] ) # just to make sure we have no zero variance
  HC = abs(F-pval)/sqrt(v)

  return( HC )
}


# determine HC decision threshold from p-values

# pval:    p-values
# alpha0:  look only at a fraction alpha0 of the p-values (default: 1)
# plot:    produce a plot of HC curve

hc.thresh = function(pval, alpha0=1, plot=TRUE)
{
  spval = sort(pval)    # sort p-values for plotting
  hcval = hc.score(spval)

  alpha0.idx = 1:ceiling(alpha0*length(hcval)) # the fraction of HC values to maximize over

  hcstat.idx = which.max( hcval[alpha0.idx] ) # idx of maximum HC 
  hcstat      = hcval[hcstat.idx]             # maximum HC
  hcstat.pval = spval[hcstat.idx]             # pval of maximum HC 

  if (plot)
  {
     plot(spval[alpha0.idx], hcval[alpha0.idx], type="l", xlab="ordered p-values", ylab="HC score")
    lines( c(hcstat.pval, hcstat.pval), c(0, hcstat)  , col=2)
  }

  return(hcstat.pval)
}


