
### slm.R  (2014-11-22)
###
###    Fit regression coefficients by plugin of (shrinkage or empirical) 
###    estimates of correlations and variances
###
### Copyright 2006-2014 Korbinian Strimmer
###
###
### This file is part of the `care' library for R and related languages.
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



slm = function(Xtrain, Ytrain, predlist, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
{
  d = dim(Xtrain)[2]

  if( missing(predlist) ) 
  {
    if(d==0)
      predlist = list(numeric(0))
    else
      predlist = list(1:d)
  }
 
  m = length(predlist)
  numpred = sapply(predlist, length)
  R2 = numeric(m)
  coeff = matrix(0, nrow=m, ncol=d+1) # plus intercept
  std.coeff = matrix(0, nrow=m, ncol=d+1) # plus intercept
  reg = matrix(0, nrow=m, ncol=2)

  modelnames = names(predlist)
  if ( is.null(modelnames) ) 
    modelnames = paste("SIZE.", numpred, sep="")
  xnames = colnames(Xtrain)
  if ( is.null(xnames) && d>0 ) 
    xnames = paste("X", 1:d, sep="")

  rownames(coeff) = modelnames  
  colnames(coeff) = c("(Intercept)", xnames)
  rownames(std.coeff) = modelnames  
  colnames(std.coeff) = c("(Intercept)", xnames)

  names(R2) = modelnames
  names(numpred) = modelnames
  colnames(reg) = c("lambda", "lambda.var")
  rownames(reg) = modelnames 
    
  for (i in 1:m)
  {
    if(verbose) cat("Determine regression coefficients for", modelnames[i], "model\n")

    # keep the best variables
    idx = predlist[[i]]

    fit = slm.internal(Xtrain[,idx, drop=FALSE], Ytrain, diagonal=diagonal, lambda=lambda, lambda.var=lambda.var, verbose=verbose)

    coeff[i, 1+idx] = fit$coefficients[-1] # predictors
    coeff[i, 1] = fit$coefficients[1] # intercept

    std.coeff[i, 1+idx] = fit$std.coefficients[-1] # predictors
    std.coeff[i, 1] = fit$std.coefficients[1] # intercept

    R2[i] = fit$R2

    reg[i,] = fit$regularization
  }

  res = list(regularization=reg, std.coefficients=std.coeff, 
             coefficients=coeff, numpred=numpred, R2=R2)
  class(res) = "slm"

  return(res)
}


make.predlist = function(ordering, numpred, name="SIZE")
{
  predlist = vector("list", length(numpred))
  names(predlist) =  paste(name, ".", numpred, sep="")
  for (i in 1:length(numpred))
    predlist[[i]] = ordering[1:numpred[i]]

  return(predlist)
}

#############

slm.internal = function(Xtrain, Ytrain, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
{
  n = dim(Xtrain)[1]
  d = dim(Xtrain)[2]
  yx = cbind(Ytrain,Xtrain)
  mu = apply(yx, 2, mean)

  # regularize the joint correlation matrix  y and x combined
  if(missing(lambda))
  {
    if (d>0)
      lambda = estimate.lambda( yx, verbose=verbose )
    else
      lambda=1
  }
  else
  {
     if(verbose) cat("Specified shrinkage intensity lambda (correlation matrix): ", round(lambda, 4), "\n")
  }
  mcor = (1-lambda)*cor(Xtrain, Ytrain) # marginal correlations

  v = var.shrink(yx, lambda.var=lambda.var, verbose=verbose)
  lambda.var = attr(v, "lambda.var")
  
  regularization = c(lambda, lambda.var)
  names(regularization) = c("lambda", "lambda.var")
  sdy = sqrt(v[1])
  sdx = sqrt(v[-1])

  if (diagonal || d==0)
  {
    bstd = mcor
  }
  else
  {
    bstd = crossprod.powcor.shrink(Xtrain, mcor, alpha=-1, lambda=lambda, verbose=FALSE)
  }

  R2 = as.vector(crossprod(mcor, bstd))   # proportion of explained variance  (diagonal=FALSE)
                                          # or sum of squared marginal correlations (diagonal=TRUE)

  b = bstd*sdy/sdx                   # regression coefficients
  a = mu[1] - sum(b * mu[-1])        # intercept

  coefficients = c(a, b)
  xnames = colnames(Xtrain)
  if ( is.null(xnames)  && d>0 ) 
    xnames = paste("X", 1:d, sep="")
  names(coefficients) = c("(Intercept)", xnames)

  std.coefficients = c(0, bstd)
  names(std.coefficients) = c("(Intercept)", xnames)

  res = list(regularization=regularization,
          std.coefficients=std.coefficients,
          coefficients=coefficients, R2=R2) 
  

  return( res )
}


