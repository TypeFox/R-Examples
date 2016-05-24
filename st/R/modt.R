### modt.R  (2015-03-21)
###
###    Moderated t Statistic
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


# Note: these function2 require the "limma" library


modt.stat = function (X, L)
{
  FUN = modt.fun(L=L)
  score = FUN(X)
  
  return( score )
}

modt.fun <- function (L)
{    
    if (missing(L)) stop("Class labels are missing!")
    L = factor(L)
    cl = levels(L)
    if (length(cl) != 2) stop("Class labels must be specified for two groups, not more or less!")

    function(X)
    {
      L = as.integer(L)
      d <- cbind(rep(1, length(L)), L)
      fit <- limma::lmFit(t(X), design=d)
      eb.out <- limma::ebayes(fit)
      modt <- -eb.out$t[,2]
  
      return(modt) 
    } 
}
