### shrinkt.R  (2013-09-01)
###
###    Shrinkage t Statistic
###
### Copyright 2006-2013 Rainer Opgen-Rhein, Verena Zuber and Korbinian Strimmer
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


shrinkt.stat = function (X, L, lambda.var, lambda.freqs, var.equal=TRUE, paired=FALSE, verbose=TRUE)
{
  if (paired)
  {
    X = pvt.pairX(X, L)
    L = rep("paired", dim(X)[1])
  }

  FUN = shrinkt.fun(L=L, lambda.var=lambda.var, lambda.freqs=lambda.freqs,
          var.equal=var.equal, verbose=verbose)

  score = FUN(X)
  
  return( score )
}

shrinkt.fun = function (L, lambda.var, lambda.freqs, var.equal=TRUE, verbose=TRUE)
{
   # shrink cat with lambda=1 equals shrink t
   return( shrinkcat.fun(L=L, lambda=1, lambda.var=lambda.var, lambda.freqs=lambda.freqs,
             var.equal=var.equal, verbose=verbose) )
}


