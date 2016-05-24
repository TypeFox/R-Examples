### pvt.svar.R  (2012-01-21)
###
###    Non-public function to compute variance shrinkage estimator 
###    
###
### Copyright 2005-12 Rainer Opgen-Rhein and Korbinian Strimmer
###
### This file is part of the `corpcor' library for R and related languages.
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


# function to compute shrinkage variance vector
#  - x: data matrix, 
#  - w:  data weights

pvt.svar = function(x, lambda.var, w, verbose)
{
  # determine correlation shrinkage intensity
  if (missing(lambda.var))
  {
    lambda.var = estimate.lambda.var(x, w, verbose)
    lambda.var.estimated=TRUE
  }
  else
  {
    if (lambda.var < 0) lambda.var = 0
    if (lambda.var > 1) lambda.var = 1
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity lambda.var (variance vector):", round(lambda.var, 4), "\n"))     
    }
    lambda.var.estimated=FALSE
  }


  # compute empirical variances 
  v = wt.moments(x, w)$var
 
  # compute shrinkage target
  target = median(v)
  
  # shrinkage estimate 
  vs = lambda.var*target + (1-lambda.var)*v
         
  attr(vs, "lambda.var") = lambda.var
  attr(vs, "lambda.var.estimated") = lambda.var.estimated
  attr(vs, "class") = "shrinkage"

  if (verbose) cat("\n")
  
  return(vs)   
}    



