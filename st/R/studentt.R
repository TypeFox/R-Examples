### studentt.R  (2013-09-01)
###
###    Student t statistic and related stuff
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


# difference of means  ("fold change")

diffmean.stat = function (X, L)
{
  FUN = diffmean.fun(L=L)
  score = FUN(X)
  
  return( score )
}

diffmean.fun = function (L)
{
    if (missing(L)) stop("Class labels are missing!")
    
    function(X)
    { 
      tmp = centroids(X, L, lambda.var=0, lambda.freqs=0, var.groups=FALSE, verbose=FALSE)
      
      # differences between the two groups
      diff = tmp$means[,1]-tmp$means[,2]
      
      return(diff)
    }
}


# student t statistic
studentt.stat = function (X, L, var.equal=TRUE, paired=FALSE)
{
  if (paired)
  {
    X = pvt.pairX(X, L)
    L = rep("paired", dim(X)[1])
  }
  
  FUN = studentt.fun(L=L, var.equal=var.equal)
  score = FUN(X)
 
  return( score )
}

studentt.fun = function (L, var.equal=TRUE)
{
   # shrinkage cat with lambda=1 and lambda.var=0 equals conventional t-score
   return( shrinkcat.fun(L=L, lambda=1, lambda.var=0, lambda.freqs=0, var.equal=var.equal, verbose=FALSE) )
}



