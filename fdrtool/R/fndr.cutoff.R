### fndr.cutoff.R  (2008-03-05)
###
###     Find Approximate Cutoff Point by an FNDR Criterion
###
### Copyright 2007-08 Korbinian Strimmer 
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


# find approximate cutoff based on false nondiscovery rate
fndr.cutoff <- function(x, statistic=c("normal", "correlation", "pvalue", "studentt"))
{
  statistic <- match.arg(statistic)
  nm = get.nullmodel(statistic)

  if(statistic=="pvalue")
  {
    ax = 1-x
    F0 <- function(zeta) return( nm$F0(zeta, sc.param) )
  }
  else
  {
    ax = abs(x)
    F0 <- function(zeta) return( 2*nm$F0(zeta, sc.param)-1 )
  }
 
  fndrfunc = function(x)
  {
    F.x = sum(ax < x)/length(ax)
    
    if (F.x == 0)
    { 
      FNDR.x = 0
    }
    else
    {
      FNDR.x = max(0, (F.x-e0.guess*F0(x)) / F.x)
    }

    return(FNDR.x)
  }


  # first, find approximate null model ("guess")
  g  = approximate.fit(x, statistic)
  sc.param = g$param
  e0.guess = g$eta0

  # second, find cutoff such that fndr is as small as possible
  MAXPCT0 = 0.99 # never use all the data
  zeta0 = as.double(quantile(ax, probs=min(MAXPCT0, e0.guess)))
  #cat("DEBUG: zeta0 =", zeta0, "\n")

  fndr2 = fndrfunc(zeta0)
  fndr1 = fndrfunc(0.9*zeta0)
  while ( fndr1 < fndr2 )
  {
    zeta0 = zeta0 *0.9
    #cat("DEBUG: zeta0 =", zeta0, "\n")
    fndr2 = fndr1
    fndr1 = fndrfunc(0.9*zeta0)
  }

  if(statistic == "pvalue") x0 = 1-zeta0
  else x0 = zeta0 
  names(x0) = "cutoff"

  rm(ax)

  return(x0)
}
