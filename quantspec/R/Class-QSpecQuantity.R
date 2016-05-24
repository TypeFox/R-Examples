#' @include generics.R
NULL

################################################################################
#' Class for a Quantile Spectral Estimator.
#'
#' \code{QSpecQuantity} is an S4 class that provides a common interface to
#' objects that are of the functional form
#' \eqn{f^{j_1, j_2}(\omega; x_1, x_2)}{f^{j1,j2}(w; x1, x2)},
#' where \eqn{j_1, j_2}{j1, j2} are indices denoting components of a time series
#' or process, \eqn{\omega}{w} is a frequency parameter and
#' \eqn{x_1, x_2}{x1, x2} are level parameters. For each combination of
#' parameters a complex number can be stored.
#' Examples for objects of this kind currently include the quantile (i. e.,
#' Laplace or copula) spectral
#' density kernel [cf. \code{\link{QuantileSD}} for an implementation], an
#' integrated version of the quantile spectral density kernels
#' [cf. \code{\link{IntegrQuantileSD}} for an implementation], and
#' estimators of it [cf. \code{\link{QuantilePG}} and \code{\link{SmoothedPG}}
#' for implementations].
#'
#' @name   QSpecQuantity-class
#' @aliases QSpecQuantity
#' @exportClass QSpecQuantity
#'
#' @keywords S4-classes
#'
#' @slot values The array holding the values
#' 							\eqn{f^{j_1, j_2}(\omega; x_1, x_2)}{f^{j1,j2}(w; x1, x2)}.
#' @slot frequencies The frequencies \eqn{\omega}{w} for which the values are
#'                    available.
#' @slot levels A list of vectors containing the levels \eqn{x_i}{xi} serving
#'               as argument for the estimator.
#'
################################################################################

setClass(
    Class = "QSpecQuantity",
    representation=representation(
        values = "array",
        frequencies = "numeric",
        levels = "list"        # a list of "numeric"s
    )
)

################################################################################
#' Get attribute \code{frequencies} from a \code{QSpecQuantity}.
#'
#' @name getFrequencies-QSpecQuantity
#' @aliases getFrequencies,QSpecQuantity-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QSpecQuantity} from which to get the \code{frequencies}.
#' @return Returns the frequencies attribute, as a vector of real numbers.
#'
#' @examples
#' qPG  <- quantilePG(rnorm(10), levels.1=c(0.25,0.5))
#' freq <- getFrequencies(qPG)
################################################################################
setMethod(f = "getFrequencies",
    signature = "QSpecQuantity",
    definition = function(object) {
      return(object@frequencies)
    }
)

################################################################################
#' Get attribute \code{levels} from a \code{QSpecQuantity}.
#'
#' If the optional parameter \code{j} is supplied, then the \code{j}th vector of
#' levels will be returned, a list with all vectors otherwise.
#'
#' @name getLevels-QSpecQuantity
#' @aliases getLevels,QSpecQuantity-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QSpecQuantity} from which to get the \code{levels}.
#' @param j Index pointing to a set of levels in the list; optional.
#'
#' @return Returns levels attribute, as a vector of real numbers.
#'
#' @examples
#' qPG         <- quantilePG(rnorm(10), levels.1=c(0.25,0.5))
#' levels.list <- getLevels(qPG)
#' levels.1    <- getLevels(qPG,1)
################################################################################
setMethod(f = "getLevels",
    signature = "QSpecQuantity",
    definition = function(object,j) {
      if (missing("j")) {
        return(object@levels)
      } else {
        if (!(j==1 | j==2)) {
          stop("Index needs to be either 1 or 2.")
        } else {
          return(object@levels[[j]])
        }
      }
    }
)

setMethod(f = "show",
    signature = "QSpecQuantity",
    definition = function(object) {

    values <- getValues(object, frequencies = object@frequencies,
        levels.1=object@levels[[1]], levels.2=object@levels[[2]])
    J <- length(object@frequencies)
    K1 <- length(object@levels[[1]])
    K2 <- length(object@levels[[2]])
    
    if (length(dim(values))==4) {
      B <- dim(values)[4]
      D <- 1
    } else { # D > 1
      B <- dim(values)[6]
      D <- dim(values)[2]
    }

    cat(paste("\n",class(object)," (J = ",J,", D = ",D,", K1 = ",K1,", K2 = ",K2,", B+1 = ",B,")\n", sep=""))

    if (J <= 7) {
      cat("Frequencies: ", object@frequencies,"\n")
    } else {
      cat("Frequencies: ", object@frequencies[1:4],"..",object@frequencies[(J-2):J],"\n")
    }

    if (K1 <= 10) {
      cat("Levels 1   : ", object@levels[[1]],"\n")
    } else {
      cat("Levels 1   : ", object@levels[[1]][1:5],"..",object@levels[[1]][(K1-4):K1],"\n")
    }

    if (K2 <= 10) {
      cat("Levels 2   : ", object@levels[[2]],"\n")
    } else {
      cat("Levels 2   : ", object@levels[[2]][1:5],"..",object@levels[[2]][(K2-4):K2],"\n")
    }

    #cat("\nValues:\n")
    if (D == 1) {
      cat("\nValues:\n")
    } else {
      cat("\nValues of first component:\n")
    }

    resultMatr <- matrix(nrow=J, ncol=K1*K2)
    cn <- rep(0,K1*K2)
    for (k1 in 1:K1) {
      for (k2 in 1:K2) {
        if (D == 1) {
          resultMatr[,k1+(k2-1)*K1] <- values[,k1,k2, 1]
        } else {
          resultMatr[,k1+(k2-1)*K1] <- values[,1,k1,1,k2, 1]
        }
        cn[k1+(k2-1)*K1] <- paste(object@levels[[1]][k1],"/",object@levels[[2]][k2], sep="")
      }
    }
    nrowShow <- min(10,nrow(resultMatr))
    ncolShow <- min(5,ncol(resultMatr))

    res <- apply(resultMatr[1:nrowShow, 1:ncolShow, drop=FALSE],c(1,2),
        function(x){complex(real=round(Re(x),3), imaginary=round(Im(x),3))})
    rownames(res) <- round(object@frequencies,3)[1:nrowShow]
    colnames(res) <- cn[1:ncolShow]

    show(res)
  }
)
