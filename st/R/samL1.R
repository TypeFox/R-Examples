### samL1.R  (2013-09-01)
###
###    Wu (2005) Improved SAM Statistic
###
### Copyright 2006-2013 Rainer Opgen-Rhein and Korbinian Strimmer
###
### This function is based in part on R code provided by Baolin Wu.
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


#### estimate Wu t-statistic (2005)

# B. Wu. 2005. Differential gene expression detection using penalized 
# linear regression model: the improved SAM statistics. 
# Bioinformatics, 21: 1565-1571

samL1.stat = function (X, L, method=c("lowess", "cor"), 
   plot=FALSE, verbose=TRUE)
{
  FUN = samL1.fun(L=L, method=method, verbose=verbose, plot=plot)
  score = FUN(X)
  
  return( score )
}

samL1.fun <- function (L, method=c("lowess", "cor"), 
   plot=FALSE, verbose=TRUE)
{
    method = match.arg(method)
   
    if (missing(L)) stop("Class labels are missing!")
     
    function(X)
    {
      tmp = centroids(X, L, lambda.var=0, lambda.freqs=0, var.groups=FALSE, verbose=verbose)
      
      # differences between the two groups
      diff = tmp$means[,1]-tmp$means[,2]
      
      # variance of diff
      n1 = tmp$samples[1]
      n2 = tmp$samples[2]
      v.diff = (1/n1 + 1/n2)*tmp$variances[,1]
      sd = sqrt(v.diff)
      
      lambda = pvt.samL1.get.lambda(diff, sd, method=method, verbose=verbose, plot=plot)
      
      # penalized t statistic
      nom = ifelse(abs(diff)>lambda, diff-sign(diff)*lambda, 0)
      den = sqrt(v.diff + lambda^2/(n1+n2-2))

      tL1 = nom/den

      return(tL1)
    }
}

## internal function

pvt.samL1.get.lambda = function(di, si, method=c("cor", "lowess"), verbose=TRUE, plot=FALSE)
{
  if (method == "lowess")
  {
    if (verbose) cat("Optimizing lambda (lowess)...\n")
  
    Lambda = seq(5, 100, length=1e2)
    rL = sapply(Lambda, function(lambda){
      mu12 = ifelse( abs(di)>lambda, di-sign(di)*lambda, 0)
      tmp = mu12/sqrt(2*si^2+lambda^2/3)
      a = lowess(si, tmp, f=2/3)$y; 
      mean((tmp[order(si)]-a)^2)/var(tmp) ## SSE/SSTO
    })
    Lam.Opt = Lambda[order(-rL)[1]] 
  
    if (plot==TRUE)
    {
      plot(Lambda, sqrt(rL), pch=20, xlab=expression(lambda), 
         ylab=expression(sqrt(SSE/SSTO)))
      abline(v=Lam.Opt, col=2, lty=2)
    }
  }
  
  if (method == "cor")
  {
    if (verbose) cat("Optimizing lambda (cor)...\n")
  
    Lambda = seq(5, 100, length=1e3)
    crL = sapply(Lambda, function(lambda){
      mu12 = ifelse( abs(di)>lambda, di-sign(di)*lambda, 0)
      tmp = mu12/sqrt(2*si^2+lambda^2/3)
      idm = order(-si)[1] ## remove the biggest si (outlier)
      return( cor(tmp[-idm],si[-idm]) )
    })
    Lam.Opt = Lambda[order(abs(crL))[1]] ## 27.14
    
    if (plot==TRUE)
    {
      plot(Lambda, abs(crL), pch=20, xlab=expression(lambda), ylab="|R|")
      abline(v=Lam.Opt, col=2, lty=2)
    }
 
  }
  
  if (verbose) cat("Estimated lamba: ", Lam.Opt, "\n")

  return(Lam.Opt)
}


