# nmfgpu4R - R binding for the nmfgpu library
# 
# Copyright (C) 2015-2016  Sven Koitka (svenkoitka@fh-dortmund.de)
# 
# This file is part of nmfgpu4R.
# 
# nmfgpu4R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# nmfgpu4R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with nmfgpu4R.  If not, see <http://www.gnu.org/licenses/>.

#' @rdname nmf
#' @export
nmf <- function(...) {
  .ensureInitialized()
  
  UseMethod("nmf")
}


#' Non-negative Matrix Factorization (NMF) on GPU 
#' 
#' Computes the non-negative matrix factorization of a data matrix \code{X} using the factorization parameter \code{r}. 
#' Multiple algorithms and initialization methods are implemented in the nmfgpu library using CUDA hardware acceleration.
#' Depending on the available hardware, these algorithms should outperform traditional CPU implementations.
#'
#' @param data Data matrix of dimension n x m with n attributes and m observations. Please note that this differs from most other data mining/machine learning algorithms!
#' 
#' @param r Factorization parameter, which affects the quality of the approximation and runtime.
#' 
#' @param algorithm Choosing the right algorithm depends on the data structure. Currently the following algorithms 
#' are implemented in the nmfgpu library:
#' \itemize{
#'  \item{\strong{mu}: Multiplicative update rules presented by Lee and Seung [2] use a purely multiplicative update
#'  and are a special form of the gradient descent algorithm (special step parameter). The implemented update rules 
#'  make use of the frobenius norm and therefore are faster than the Kullback-Leibler ones.}
#'  
#'  \item{\strong{gdcls}: Gradient Descent Constrained Least Squares (GDCLS) presented by Pauca et al [3,4] is a hybrid 
#'  algorithm. It uses a least squares solver for updating matrix \code{H} and the multiplicative update rule for \code{W}
#'  as defined by Lee and Seung [2]. Additionaly the GDCLS algorithm uses the parameter \code{lambda} from the 
#'  \code{parameters} argument to control the sparsity of the matrix \code{H}. As from the authors presented, the 
#'  sparsity parameter \code{lambda} should be between \code{0.1} and \code{0.0001}, but at least positive.}
#'  
#'  \item{\strong{als}: Alternating Least Squares (ALS) originally presented by Paatero and Tapper [1,4] uses a 
#'  least squares solver for updating both matrices \code{W} and \code{H}.}
#'  
#'  \item{\strong{acls}: Alternating Constrained Least Squares (ACLS) presented by Langville et al [4] enhances the normal 
#'  ALS algorithm by introducing sparsity parameters. Both \code{lambdaW} and \code{lambdaH} must be provided in the
#'  \code{parameters} argument and must be in the range of \code{0} and positive infinity.}
#'  
#'  \item{\strong{ahcls}: Alternating Hoyer Constrained Least Squares (AHCLS) presented by Langville et al [4] enhances 
#'  the ACLS algorithm by introducing a second set of sparsity parameters. Additionaly to \code{lambdaW} and \code{lambdaH}
#'  the sparsity parameters \code{alphaH} and \code{alphaW} must be provided in the \code{parameters} argument. Both should
#'  be set in the range of \code{0.0} and \code{1.0}, representing a percentage sparsity for each matrix. As recommended by
#'  the authors all four parameters should be set to \code{0.5} as starting values.)}
#'  \item{\strong{nsnmf}: Non-smooth Non-negative Matrix Factorization (nsNMF) presented by Pascual-Montano et al [6] is an
#'  enhancement to the multiplicative update rules [2]. With an extra parameter \code{theta} the user has control over the 
#'  influence of the model. The value should be in the range of \code{0.0} and \code{1.0} to work like intended.}
#' }
#' 
#' @param initMethod All initialization methods depend on the selected algorithm. Using the fact that a least squares type 
#' algorithm computes the matrix \code{H} in the first step, does make an initialization for \code{H} unnecessary. Therefore 
#' only the initialization method for matrix \code{W} will be executed for any least squares type algorithm.
#' \itemize{
#'  \item{\strong{CopyExisting}: Initializes the factorization matrices \code{W} and \code{H} with existing values, 
#'  which requires \code{W} and \code{H} to be set in the \code{parameters} argument. On the one hand this enables the
#'  user to chain different algorithms, for example using a fast converging algorithm for a base approximation and and a
#'  slow algorithm with better convergence properties to finish the optimization process. On the other hand the user can
#'  supply matrix intializations, which are not supported by this interface.
#'  \emph{Note}: Both \code{W} and \code{H} must have the same dimension as they would have from the passed arguments 
#'  \code{X} and \code{r}.}
#'  
#'  \item{\strong{AllRandomValues}: Initializes the factorization matrices \code{W} and \code{H} with uniformly distributed 
#'  values between \code{0.0} and \code{1.0}, where \code{0.0} is excluded and \code{1.0} is included.}
#'  
#'  \item{\strong{MeanColumns}: Initializes the factorization matrix \code{W} by computing the mean of five random data 
#'  matrix columns. The matrix \code{H} will be initialized as it would when using \code{AllRandomValues}.}
#'  
#'  \item{\strong{k-Means/Random}: Initializes the factorization matrix \code{W} by computing the k-Means cluster centers of
#'  the data matrix. The matrix \code{H} will be initialized as it would when using \code{AllRandomValues}. 
#'  This method was presented by Gong et al [5] as initialization strategy H2.}
#'  
#'  \item{\strong{k-Means/NonNegativeWTV}: Initializes the factorization matrix \code{W} by computing the k-Means cluster centers of
#'  the data matrix. The matrix \code{H} will be initialized with the product \code{t(W) \%*\% V}, but all negative values are clamped to zero. 
#'  This method was presented by Gong et al [5] as initialization strategy H4.}
#'  
#'  \item{\strong{EIn-NMF}: Initializes the factorization matrix \code{W} by computing the k-Means cluster centers of the data matrix.
#'  The matrix \code{H} will be initialized with a prefix sum equation to build weighted encoding vectors. This method was presented by Gong
#'  et al [5] as initialization strategy H5.}
#' }
#' 
#' @param seed The \code{seed} is used to initialize the random number generators for initializing the factorization matrices. 
#' Setting this argument to a fixed value ensures the same initialization of the matrices. This can be handy for benchmarking
#' different algorithms with the same initialization. 
#' 
#' @param threshold \strong{First convergence criterion:} The \code{threshold} is used to determine if the algorithm has 
#' converged to a local minimum by checking the difference between the last frobenius norm and the current one. If it is
#' below the \code{threshold}, then the algorithm is assumed to be converged.
#' 
#' @param maxiter \strong{Second convergence criterion:} If the first convergence criterion is not reached, but a maximum
#' number of iterations, the execution of the algorithm will be aborted.
#' 
#' @param runs Performs the specified amount of runs and stores the best factorization result.
#' 
#' @param parameters A list of additional parameters, which are required by some \code{algorithm} and \code{initMethod}
#' values.
#' 
#' @param useSinglePrecision By default R only knows about double precision numerical data types. If this parameter is set to true,
#' then the algorithm will convert the double precision data to single precision for computation. The result will be converted back to 
#' double precision data.
#' 
#' @param verbose By default information about the factorization process and current status will be written to the console. For silent execution \code{verbose=T} may be passed, preventing any output besides error messages.
#' 
#' @param ssnmf Internal flag (Don't use it)
#' 
#' @param ... Other arguments
#' 
#' @return If the factorization process was successful, then a list of the following values will be returned otherwise NULL:
#' \tabular{ll}{
#'  \code{W} \tab Factorized matrix W with n attributes and r basis features of the data matrix.\cr
#'  \code{H} \tab Factorized matrix H with r mixing vectors for m data entries in the data matrix.\cr
#'  \code{Frobenius} \tab Contains the frobenius norm of the factorization at the end of algorithm execution.\cr
#'  \code{RMSD} \tab Contains the root-mean-square deviation (RMSD) of the factorization at the end of algorithm execution.\cr
#'  \code{ElapsedTime} \tab Contains the elapsed time for initialization and algorithm execution.\cr
#'  \code{NumIterations} \tab Number of iterations until the algorithm had converged.\cr
#' }
#' 
#' @examples
#' \dontrun{
#' # Initialize the library
#' library(nmfgpu4R)
#' nmfgpu4R.init()
#' 
#' # Create dummy data
#' data <- runif(256*1024)
#' dim(data) <- c(256, 1024)
#' 
#' # Compute several factorization models
#' result <- nmf(data, 128, algorithm="mu", initMethod="K-Means/Random", maxiter=500)
#' result <- nmf(data, 128, algorithm="mu", initMethod="CopyExisting", 
#'                  parameters=list(W=result$W, H=result$H), maxiter=500)
#' result <- nmf(data, 128, algorithm="gdcls", maxiter=500, parameters=list(lambda=0.1))
#' result <- nmf(data, 128, algorithm="als", maxiter=500)
#' result <- nmf(data, 128, algorithm="acls", maxiter=500, 
#'                  parameters=list(lambdaH=0.1, lambdaW=0.1))
#' result <- nmf(data, 128, algorithm="ahcls", maxiter=500, 
#'                  parameters=list(lambdaH=0.1, lambdaW=0.1, alphaH=0.5, alphaW=0.5))
#' result <- nmf(data, 128, algorithm="nsnmf", maxiter=500, parameters=list(theta=0.25))
#' 
#' # Compute encoding matrices for training and test data
#' set.seed(42)
#' idx <- sample(1:nrow(iris), 100, replace=F)
#' data.train <- iris[idx,]
#' data.test <- iris[-idx,]
#' 
#' model.nmf <- nmf(t(data.train[,-5]), 2)
#' encoding.train <- t(predict(model.nmf))
#' encoding.test <- t(predict(model.nmf, t(data.test[,-5])))
#' 
#' plot(encoding.train, col=data.train[,5], pch=1)
#' points(encoding.test, col=data.test[,5], pch=4)
#' }
#' 
#' @references
#' \enumerate{
#'  \item{P. Paatero and U. Tapper, "Positive matrix factorization: A non-negative factor model with optimal utilization of error estimates of data values", Environmetrics, vol. 5, no. 2, pp. 111-126, 1994.}
#'  \item{D. D. Lee and H. S. Seung, "Algorithms for non-negative matrix factorization", in Advances in Neural Information Processing Systems 13 (T. Leen, T. Dietterich, and V. Tresp, eds.), pp. 556-562, MIT Press, 2001.}
#'  \item{V. P. Pauca, J. Piper, and R. J. Plemmons, "Nonnegative matrix factorization for spectral data analysis", Linear Algebra and its Applications, vol. 416, no. 1, pp. 29-47, 2006. Special Issue devoted to the Haifa 2005 conference on matrix theory.}
#'  \item{A. N. Langville, C. D. Meyer, R. Albright, J. Cox, and D. Duling, "Algorithms, initializations, and convergence for the nonnegative matrix factorization", CoRR, vol. abs/1407.7299, 2014.}
#'  \item{L. Gong and A. Nandi, "An enhanced initialization method for non-negative matrix factorization", in 2013 IEEE International Workshop on Machine Learning for Signal Processing (MLSP), pp. 1-6, Sept 2013.}
#'  \item{A. Pascual-Montano, J. M. Carazo, K. Kochi, D. Lehmann and R. D. Pascual-Marqui "Nonsmooth nonnegative matrix factorization (nsNMF)", in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 28, pp. 403-415, 2006}
#' }
#' 
#' @rdname nmf
#' @method nmf default
#' @export
nmf.default <- function(data, r, algorithm="mu", initMethod="AllRandomValues", seed=floor(runif(1, 0, .Machine$integer.max)), threshold=0.1, maxiter=2000, runs=1, parameters=NULL, useSinglePrecision=F, verbose=T, ssnmf=F, ...) {
  if(r > dim(data)[2] && !ssnmf) {
    stop("Factorization parameter r must be less or equal the number of columns of the data matrix!")
  }
  
  # If it is a data.frame then try to convert it to a numeric matrix
  if(is.data.frame(data)) {
    if(all(sapply(data, is.numeric))) {
      data <- as.matrix(data)
    } else {
      stop("Data frame does contain non-numeric variables")
    }
  }
  
  validateMatrix(data)
  
  if(isS4(data)) {
    if(is.matrix.SparseM(data) || is.matrix.Matrix(data)) {
      if(useSinglePrecision) {
        result <- adapterComputeSinglePrecisionSparse(algorithm, initMethod, data, r, seed, threshold, maxiter, runs, parameters, verbose, ssnmf)
      } else {
        result <- adapterComputeDoublePrecisionSparse(algorithm, initMethod, data, r, seed, threshold, maxiter, runs, parameters, verbose, ssnmf)        
      }
    } else {
      stop("Unknown S4 class as data matrix!")
    }
  } else {
    if(useSinglePrecision) {
      result <- adapterComputeSinglePrecision(algorithm, initMethod, data, r, seed, threshold, maxiter, runs, parameters, verbose, ssnmf)
    } else {
      result <- adapterComputeDoublePrecision(algorithm, initMethod, data, r, seed, threshold, maxiter, runs, parameters, verbose, ssnmf)
    }
  }
  
  if(!is.null(result)) {
    # Copy row and column names from data matrix and generate new names for basis vectors
    rownames(result$W) <- rownames(data)
    colnames(result$W) <- paste("r", 1:r, sep="")
    rownames(result$H) <- colnames(result$W)
    colnames(result$H) <- colnames(data)
    
    # Set class
    class(result) <- "nmfgpu"
    
    # Prepare for SSNMF
    result$call <- list()
    result$call$arglist <- list(...)
    result$call$arglist$r <- r
    result$call$data.nrow <- nrow(data)
  }
  
  return(result)
}

#' @param formula Formula object with no intercept and labels for selected attributes. Note that die labels are 
#' selected from the rows instead of the columns, because NMF expects rows to be attributes.
#'
#' @rdname nmf
#' @method nmf formula
#' @export
nmf.formula <- function(formula, data, ...) { 
  labels <- validateFormulaAndGetLabels(formula, data)
  
  result <- nmf.default(data[labels,], ...)
  if(!is.null(result)) {
    result$call$labels <- labels
  }
  return(result)
}

# @rdname nmf
# @method residuals nmfgpu
# @export
#residuals.nmfgpu <- function(object, ...) {
#  return(object$Frobenius)
#}

#' @rdname nmf
#' @method fitted nmfgpu
#' @export
fitted.nmfgpu <- function(object, ...) {
  return(object$W %*% object$H)
}

#' @param object Object of class "\code{nmfgpu}"
#' @param newdata New data matrix compatible to the training data matrix, for computing the corresponding mixing matrix.
#' 
#' @rdname nmf
#' @method predict nmfgpu
#' @export
predict.nmfgpu <- function(object, newdata, ...) {
  .ensureInitialized()
  
  # If no new data is provided then return the encoding matrix of the training data
  if(missing(newdata)) {
    return(object$H)
  }
  
  # If new data is a data.frame, then try to convert it
  if(is.data.frame(newdata)) {
    newdata <- as.matrix(newdata)
  }
  
  # If a formula was provided to the training data, then restrict newdata to the same features
  if("labels" %in% names(object$call)) {
    if(!all(object$call$labels %in% rownames(newdata))) {
      stop("Data matrix does not contain the same features as the training data!")
    } else {
      newdata <- newdata[,object$call$labels]
    }
  }
  
  # Check if new data has the same amount of features as the training data
  if(is.matrix(newdata) || is.matrix.SparseM(newdata) || is.matrix.Matrix(newdata)) {
    if(nrow(newdata) != object$call$data.nrow) {
      stop("New data matrix does not have the same amount of features as the training data!")
    }
  }
  
  # Adjust original arguments to support SSNMF
  args <- object$call$arglist
  args$initMethod <- "AllRandomValues"
  if("parameters" %in% names(args)) {
    args$parameters$W <- object$W
  } else {
    args$parameters <- list(W=object$W)
  }
  
  # Compute the NMF for the new data using the basis vectors from the training data
  args$data <- newdata
  args$runs <- 1 # Encoding with fixed basis vectors is a convex optimization!
  args$ssnmf <- T
  
  result <- do.call(nmf.default, args)
  
  if(!is.null(result)) {
    return(result$H)
  } else {
    return(NULL)
  }
}