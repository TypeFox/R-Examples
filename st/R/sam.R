### sam.R  (2015-03-21)
###
###    SAM t Statistic
###
### Copyright 2006-2015 Rainer Opgen-Rhein and Korbinian Strimmer
###
###
### This file is part of the `st' library for R and related languages.
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


# Note: these functions require the "samr" library

sam.stat = function (X, L)
{
  FUN = sam.fun(L=L)
  score = FUN(X)
  
  return( score )
}

sam.fun <- function(L)
{
    if (missing(L)) stop("Class labels are missing!")
    L = factor(L)
    cl = levels(L)
    if (length(cl) != 2) stop("Class labels must be specified for two groups, not more or less!")
  
   function(X)
    {
      dd = list(x=t(X),y=as.integer(L), logged2=TRUE)
      out = samr::samr(dd, resp.type="Two class unpaired", nperms=1)
  
      return(out$tt) # SAM test statistic
    } 
}
