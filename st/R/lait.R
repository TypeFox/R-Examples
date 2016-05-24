### lait.R  (2013-09-01)
###
###    Correlation-Predicted t-Statistic of Lai
###
### Copyright 2009-13 Verena Zuber and Korbinian Strimmer
###
### Parts of the code are adopted from example code by Yinglai Lai
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


lait.stat = function (X, L, f=0.2, verbose=TRUE)
{
  FUN = lait.fun(L=L, f=f, verbose=verbose)
  score = FUN(X)
  
  return( score )
}

lait.fun = function (L, f=0.2, verbose=TRUE)
{
  if (missing(L)) stop("Class labels are missing!")
 
  function(X)
  { 
    tmp = centroids(X, L, lambda.var=0, lambda.freqs=0, var.groups=FALSE, 
            centered.data=TRUE, verbose=verbose)
      
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

    if (verbose) cat("\nComputing Lai's correlation-predicted t-statistics\n")
    m = length(t)
    score = sapply(1:m, lai.tscore, t, R, f=f)
      
    return( score )
  }
}


lai.tscore = function(gene, tscore, corr, f=0.2, plot=FALSE)
{
  cscore = corr[gene,]
  s = sign(cscore)
  x = s*cscore
  y = s*c(tscore)
  
  mod = lowess(x[-gene], y[-gene], f)
  m = length(tscore)
  coe = lsfit(x=mod$x[(m-2):(m-1)], y=mod$y[(m-2):(m-1)])$coe
  ynew = coe[1] + coe[2]
 
  if (plot)
  {
    plot(x, y, xlab="absolute correlation", ylab="t statistic")
    lines(mod, col=2)
    points(1, ynew, col=2) 
  }

  return(as.double(ynew))
}



