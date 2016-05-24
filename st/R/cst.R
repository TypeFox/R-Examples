### cst.R  (2013-09-01)
###
###    Correlation-Shared t-Statistic of Tibshirani-Wassermann
###
### Copyright 2009-2013 Korbinian Strimmer
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


cst.stat = function (X, L, verbose=TRUE)
{
  FUN = cst.fun(L=L, verbose=verbose)
  score = FUN(X)
  
  return( score )
}

cst.fun = function (L, verbose=TRUE)
{
    if (missing(L)) stop("Class labels are missing!")
 
    function(X)
    { 
      tmp = centroids(X, L, var.groups=FALSE, centered.data=TRUE,
              lambda.var=0, lambda.freqs=0, verbose=verbose)
      
      # differences between the two groups
      diff = tmp$means[,1]-tmp$means[,2]
      
      # standard error of diff
      n1 = tmp$samples[1]
      n2 = tmp$samples[2]
      v =  tmp$variances[,1]  # pooled variance
      sd = sqrt( (1/n1 + 1/n2)*v )
      
      # pooled empirical correlation matrix
      R = cor(tmp$centered.data) 	

      # t statistic
      t = diff/sd

      if (verbose) cat("\nComputing correlation-shared t-statistics\n")
 
      # compute correlation-shared t-statistic
      p = length(t)
      cst.vec = numeric(p)
      for (i in 1:p)  
      {
        idx = order(R[i,], decreasing=TRUE)         # sort other genes by correlation
        nonneg = sum(R[i,] >= 0)                    # number of nonnegatively associated genes
        y = cumsum(abs(t[idx[1:nonneg]]))/1:nonneg  # find maximum average 
        cst.vec[i] = max(y)*sign(t[i])
      }

      return(cst.vec)
    }
}

