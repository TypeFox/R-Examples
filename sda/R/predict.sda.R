### predict.sda.R  (2013-11-21)
###
###    Shrinkage discriminant analysis (prediction)
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


predict.sda = function (object, Xtest, verbose = TRUE, ...) 
{
    if (missing(object)) {
        stop("A sda fit object must be supplied.")
    }
    if (missing(Xtest)) {
        stop("A new data to predict must be supplied.")
    }
    if (!is.matrix(Xtest)) 
        stop("Test data must be given as matrix!")
    ntest = nrow(Xtest)

    alpha = object$alpha
    cl.count = length(alpha)

    if (ncol(Xtest) != ncol(object$beta)) 
        stop("Different number of predictors in sda object (", 
            ncol(object$beta), ") and in test data (", ncol(Xtest), 
            ")", sep = "")
     
    beta = object$beta
    if (verbose) 
        cat("Prediction uses", ncol(beta), "features.\n")
      
    probs = t(tcrossprod(beta, Xtest) + alpha)
    probs = exp(probs - max.col.value(probs))  #probs = exp(probs - apply(probs, 1, max))
    probs = zapsmall( probs / rowSums(probs) )

    yhat = max.col(probs) # yhat = apply(probs, 1, which.max)

    attr(yhat, "levels") = names(alpha)
    class(yhat) = "factor"
    colnames(probs) = names(alpha)
    rownames(probs) = rownames(Xtest)
    return(list(class = yhat, posterior = probs))
}

# by Sebastian Gibb
max.col.value = function(x)
{
  return(x[cbind(1:nrow(x), max.col(x, ties.method="first"))])
}



