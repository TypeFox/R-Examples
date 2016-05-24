### shrink.intensity.R  (2012-09-02)
###
###   Functions for computing the shrinkage intensity
###    
###
### Copyright 2005-2012 Juliane Sch\"afer, Rainer Opgen-Rhein, 
###                     Miika Ahdesm\"aki and Korbinian Strimmer
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


# estimate shrinkage intensity lambda.var (variance vector)
#
# input:  data matrix
#         weights of each data point
#
estimate.lambda.var = function(x, w, verbose=TRUE)
{
  n = nrow(x) 
  w = pvt.check.w(w, n)
  
  # bias correction factors
  w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xc)[1]
  h1 = 1/(1-w2)           # for w=1/n this equals the usual h1=n/(n-1)
  h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)
   
  # center input matrix
  xc = wt.scale(x, w, center=TRUE, scale=FALSE) 

  # compute empirical variances 
  #v = wt.moments(x, w)$var
  v = h1*(colSums(w*xc^2))

  # compute shrinkage target
  target = median(v)


  if (verbose) cat("Estimating optimal shrinkage intensity lambda.var (variance vector): ")

  zz = xc^2
  q1 = colSums( sweep(zz, MARGIN=1, STATS=w, FUN="*") )
  q2 = colSums( sweep(zz^2, MARGIN=1, STATS=w, FUN="*") ) - q1^2   
  numerator = sum( q2 )
  denominator = sum( (q1 - target/h1)^2 )

  if(denominator == 0) 
    lambda.var = 1
  else
    lambda.var = min(1, numerator/denominator * h1w2)
 
  if (verbose) cat(paste(round(lambda.var, 4), "\n")) 
  
  return (lambda.var)
}


# estimate shrinkage intensity lambda (correlation matrix)
#
# input:  data matrix
#         weights of each data point
#
#
# note: the fast algorithm in this function is due to Miika Ahdesm\"aki
#
estimate.lambda = function(x, w, verbose=TRUE)
{
   n = nrow(x)
   p = ncol(x)

   if (p == 1) return (1) 

   w = pvt.check.w(w, n)
   xs = wt.scale(x, w, center=TRUE, scale=TRUE) # standardize data matrix

  if (verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix): ")

   # bias correction factors
   w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xs)[1]
   h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)

   sw = sqrt(w)
   
   # direct slow algorithm
   #  E2R = (crossprod(sweep(xs, MARGIN=1, STATS=sw, FUN="*")))^2
   #  ER2 = crossprod(sweep(xs^2, MARGIN=1, STATS=sw, FUN="*"))
   #  ## offdiagonal sums
   #  sE2R = sum(E2R)-sum(diag(E2R))
   #  sER2 = sum(ER2)-sum(diag(ER2))
 
   # Here's how to compute off-diagonal sums much more efficiently for n << p
   # this algorithm is due to Miika Ahdesm\"aki
   xsw = sweep(xs, MARGIN=1, STATS=sw, FUN="*")
   xswsvd = fast.svd(xsw)
   sE2R = sum(xsw*(sweep(xswsvd$u,2,xswsvd$d^3,'*')%*%t(xswsvd$v))) - sum(colSums(xsw^2)^2) 
   remove(xsw,xswsvd) # free memory 
   xs2w = sweep(xs^2, MARGIN=1, STATS=sw, FUN="*")
   sER2 = 2*sum(xs2w[,(p-1):1] * t(apply(xs2w[,p:2, drop=FALSE],1,cumsum)))
   remove(xs2w) # free memory 
 
   #######

   denominator = sE2R
   numerator = sER2 - sE2R

   if(denominator == 0)
     lambda = 1
   else
     lambda = min(1, numerator/denominator * h1w2)

   if (verbose) cat(paste(round(lambda, 4), "\n"))

   return (lambda)
}

