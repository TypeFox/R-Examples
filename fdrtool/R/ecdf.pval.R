### ecdf.pval.R  (2007-06-15)
###
###    Estimate Empirical Density of p-Values 
###    
###
### Copyright 2007 Korbinian Strimmer
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


# empirical cumulative distribution of p-values,
# constrained such that the known fraction eta0 of null p-values 
# is taken into account
ecdf.pval <- function (x, eta0=1) 
{
    # compute empirical CDF as usual
    x = sort(x)
    n = length(x)
    if (n < 1) 
        stop("'x' must have 1 or more non-missing values")
    vals = sort(unique(x))
    F.raw = cumsum(tabulate(match(x, vals)))/n
    
    # control upper bound of F:
    # make sure that the maximum slope of (Grenander) F is eta0
    F.raw = pmin(F.raw, 1-eta0*(1-vals) ) 

    # control lower bound of F: 
    # make sure that (Grenander F) >= eta0*vals
    F.raw = pmax(F.raw, eta0*vals) 

    # if necessary add an atom at 1 to make it a proper CDF
    if (vals[length(vals)] != 1)
    {
       F.raw = c(F.raw, 1)
       vals = c(vals, 1)
    }

    # if necessary also add an atom at 0 with weight zero to get support [0,1]
    if (vals[1] != 0)
    {
       F.raw = c(0, F.raw)
       vals = c(0, vals)
    }

    # finally, modify F such that the last slope of the Grenander F 
    # is *exactly* eta0
    i = length(vals)-1
    F.raw[i] = 1-eta0*(1-vals[i])
    
    rval <- approxfun(vals, F.raw, 
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) = c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    rval
}


