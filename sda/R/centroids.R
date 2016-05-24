### centroids.R  (2013-08-31)
###
###    Group centroids and (pooled) variances
###
### Copyright 2008-2013 Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
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



centroids = function(x, L, lambda.var, lambda.freqs,
  var.groups=FALSE, centered.data=FALSE, verbose=TRUE)
{
  if (!is.matrix(x)) stop("Input x must be a matrix!")
  p = ncol(x)
  n = nrow(x)

  if (length(L) != n) stop("For each sample there must be a class label!")   
  
  groups = pvt.groups(L)
  samples = groups$samples
  cl.count = groups$cl.count
  cl.names = groups$cl.names

  if( missing(lambda.var) )
    auto.shrink = TRUE
  else
  {
    auto.shrink = FALSE

    if(var.groups)
    {
      if (length(lambda.var) == 1)
        specified.lambda.var = rep(lambda.var[1], (cl.count+1))
      else
      {
        if (length(lambda.var) != (cl.count+1) )
          stop("You need to specify a vector with ", (cl.count+1), " shrinkage intensities, one for each group variance plus one for the pooled variance!")
        specified.lambda.var = lambda.var
      }
    }
    else
      specified.lambda.var = lambda.var[1]

    ifelse(specified.lambda.var > 1, 1, specified.lambda.var)
    ifelse(specified.lambda.var < 0, 0, specified.lambda.var)

  }


  if (verbose)
  {
    cat("Number of variables:", p, "\n")
    cat("Number of observations:", n, "\n")
    cat("Number of classes:", cl.count, "\n\n")
  }

  # estimate class frequencies
  lambda.freqs.estimated = missing(lambda.freqs)
  freqs = freqs.shrink( samples, lambda.freqs=lambda.freqs, verbose=verbose )
  attr(freqs, "lambda.freqs.estimated") = lambda.freqs.estimated
 
  
  # setup arrays
  mu = array(0, dim=c(p, cl.count+1))
  colnames(mu) = c(cl.names, "(pooled)")
  rownames(mu) = colnames(x)

  xc = array(0, dim=c(n,p))  # storage for centered data
  rownames(xc) = rownames(x)
  colnames(xc) = colnames(x)
 
  if(var.groups)
  {
    v = array(0, dim=c(p, cl.count+1)) # storage for variances
    attr(v, "lambda.var") = numeric(cl.count+1)
    colnames(v) = c(cl.names, "(pooled)")
    rownames(v) = colnames(x)
  }
  else
  {
    v = array(0, dim=c(p, 1)) # store only pooled variances
    attr(v, "lambda.var") = numeric(1)
    colnames(v) = c("(pooled)")
    rownames(v) = colnames(x)
  }  

  # compute means and variance in each group
  mu.pooled = rep(0, p)
  for (k in 1:cl.count)
  {
     idx = groups$idx[,k]
     Xk = x[ idx, ,drop = FALSE]
     mu[,k] = colMeans(Xk)
     mu.pooled = mu.pooled + freqs[k]*mu[,k]

     xc[idx,] = sweep(Xk, 2, mu[,k]) # center data

     if (var.groups)
     {
        if(verbose) cat("Estimating variances (class #", k, ")\n", sep="")
        if (auto.shrink)
        {
          vs = var.shrink(Xk, verbose=verbose)
        }
        else
        {
          vs = var.shrink(Xk, lambda.var=specified.lambda.var[k], verbose=verbose)
        }
        v[,k] = as.vector(vs)
        attr(v, "lambda.var")[k] = attr(vs, "lambda.var")  
     }
  }
  mu[,cl.count+1] = mu.pooled
 
  
  # compute variance
  if (verbose) cat("Estimating variances (pooled across classes)\n")
  if (var.groups)
  {
    if (auto.shrink)
      v.pool = var.shrink(xc, verbose=verbose)
    else
      v.pool = var.shrink(xc, lambda.var=specified.lambda.var[cl.count+1], verbose=verbose)

    v[,cl.count+1] = v.pool*(n-1)/(n-cl.count) # correction factor
    attr(v, "lambda.var")[cl.count+1] = attr(v.pool, "lambda.var")
  }  
  else
  {
    if (auto.shrink)
      v.pool = var.shrink(xc, verbose=verbose)
    else
      v.pool = var.shrink(xc, lambda.var=specified.lambda.var[1], verbose=verbose)
    v[,1] = v.pool*(n-1)/(n-cl.count) # correction factor
    attr(v, "lambda.var")[1] = attr(v.pool, "lambda.var")
  }

  if (auto.shrink)
    attr(v, "lambda.var.estimated") = TRUE
  else
    attr(v, "lambda.var.estimated") = FALSE

  if(centered.data == FALSE) xc=NULL
  
  ##

  return( list(samples=samples, freqs= freqs, means=mu, variances=v, 
     centered.data=xc))
}



## private function ##

pvt.groups = function(L)
{
   y = factor(L)  # note that this creates new levels (in contrast to as.factor)
     
   cl.names = levels(y)
   cl.count = length(cl.names)

   idx = array(FALSE, dim = c(length(y), cl.count))
   colnames(idx) = cl.names
   nn = integer(cl.count)
   names(nn) = cl.names
   
   for (k in 1:cl.count)
   {
      idx[,k] = ( y == cl.names[k] )
      nn[k] = sum(idx[,k])
   }

   return( list(idx=idx, samples=nn, cl.count=cl.count, cl.names=cl.names) )
}

