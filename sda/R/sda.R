### sda.R  (2013-11-19)
###
###    Shrinkage discriminant analysis (training the classifier)
###
### Copyright 2008-13 Miika Ahdesmaki and Korbinian Strimmer
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


sda = function(Xtrain, L, lambda, lambda.var, lambda.freqs, diagonal=FALSE, verbose=TRUE)
{
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(L)) stop("Class labels are missing!")

  # shrinkage intensities
  regularization = rep(NA, 3)
  names(regularization) = c("lambda", "lambda.var", "lambda.freqs")
  regularization[1] = 1 # for diagonal=TRUE


  tmp = centroids(Xtrain, L, lambda.var, lambda.freqs, var.groups=FALSE, centered.data=TRUE, verbose=verbose)
  
  cl.count = ncol(tmp$means)-1    # number of classes
  n = sum(tmp$samples)            # number of samples
  p = nrow(tmp$means)             # number of features

  mu = tmp$means[,1:cl.count, drop=FALSE]
  mup = tmp$means[,cl.count+1, drop=FALSE]
  sc = sqrt(tmp$variances[,1])  
  regularization[2] = attr(tmp$variances, "lambda.var")[1]

  # class frequencies
  freqs = tmp$freqs
  regularization[3] = attr(freqs, "lambda.freqs")
  attr(freqs, "lambda.freqs") = NULL
  attr(freqs, "lambda.freqs.estimated") = NULL
 
  xc = tmp$centered.data # to compute inverse pooled correlation matrix

  rm(tmp)

  
  ############################################################# 
  # compute coefficients for prediction 
  #############################################################
  
  # prediction weights
  pw = array(0, dim=c(p, cl.count) )
  colnames(pw) = colnames(mu)
  rownames(pw) = rownames(mu)

  for (k in 1:cl.count)
  {
    diff = mu[,k]-mup  
    pw[,k] = diff/sc
  }

  if(!diagonal)
  {
    if(verbose) cat("\nComputing inverse correlation matrix (pooled across classes)\n")
    pw = crossprod.powcor.shrink(xc, pw, alpha=-1, lambda=lambda, verbose=FALSE)
    regularization[1] = attr(pw, "lambda")
    lambda.estimated = attr(pw, "lambda.estimated")
    attr(pw, "lambda") = NULL
    attr(pw, "lambda.estimated") = NULL
    if(verbose)
      if(lambda.estimated)  
        cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
                    round(regularization[1], 4), "\n")
      else
         cat("Specified shrinkage intensity lambda (correlation matrix):", 
                    round(regularization[1], 4), "\n")
  }

  for (k in 1:cl.count)
  {
    pw[,k] = pw[,k]/sc
  }

  alpha = log(freqs)
  for (k in 1:cl.count) {
    refk = (mu[,k]+mup)/2
    alpha[k] = alpha[k]-crossprod(pw[,k], refk) 
  }


  ############################################################# 

  out = list(regularization=regularization, freqs=freqs, 
             alpha=alpha, beta=t(pw))
  class(out)="sda"

  return (out)
}

