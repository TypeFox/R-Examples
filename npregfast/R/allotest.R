#' Bootstrap based test for testing an allometric model
#'@description Bootstrap-based procedure that tests whether the data 
#' can be modelled by an allometric model.
#'@param formula An object of class \code{formula}: a sympbolic description
#' of the model to be fitted.
#'@param data A data frame argumment or matrix containing the model response variable
#' and covariates required by the \code{formula}.
#'@param nboot Number of bootstrap repeats.
#'@param kbin Number of binning nodes over which the function
#' is to be estimated.
#'@param seed Seed to be used in the bootstrap procedure.
#'@details In order to facilitate the choice of a model appropriate
#' to the data while at the same time endeavouring to minimise the 
#' loss of information,  a bootstrap-based procedure, that test whether the 
#' data can be modelled by an allometric model, was developed.  Therefore,
#' \code{allotest} tests the null hypothesis of an allometric model taking 
#' into account the logarithm of the original variable
#'  (\eqn{X^* = log(X)}{} and \eqn{Y^* =log (Y)}{}). 
#' Based on a general model of the type 
#' \deqn{Y^*=m(X^*)+\varepsilon}{} 
#' the aim here is to test the null hypothesis of an allometric model 
#' \deqn{H_0 = m(x^*) =  a^*+ b^* x^*}{} 
#' \eqn{vs.}{} the general hypothesis 
#' \eqn{H_1}{}, with \eqn{m}{} 
#' being an unknown nonparametric function; or analogously,
#' \deqn{H_1: m(x^*)= a^*+ b^* x^* + g(x^*)}{}
#' with \eqn{g(x^*)}{} being an unknown function not equal to zero. 
#' 
#' To implement this test we have used the wild bootstrap.
#' 
#' 
#' 
#'@return An object is returned with the following elements:
#' \item{statistic}{the value of the test statistic.}
#' \item{value}{the p-value of the test.}
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@references 
#' Sestelo, M. and Roca-Pardinas, J. (2011). A new approach to estimation of 
#' length-weight relationship of \eqn{Pollicipes}  \eqn{pollicipes} (Gmelin, 1789)
#' on the Atlantic coast of Galicia (Northwest Spain): some aspects of its 
#' biology and management. Journal of Shellfish Research, 30 (3), 939--948.
#' 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#' 
#'@examples
#' library(npregfast)
#' data(barnacle)
#' allotest(DW ~ RC, data = barnacle, nboot = 100, seed = 130853)
#' 
#' @useDynLib npregfast allotest_
#' @importFrom stats na.omit runif
#' @export





allotest <- function(formula, data = data, nboot = 500, kbin = 200, 
                     seed = NULL) {
  
  ffr <- interpret.frfastformula(formula, method = "frfast")
  varnames <- ffr$II[2, ]
  aux <- unlist(strsplit(varnames, split = ":"))
  varnames <- aux[1]
  namef <- aux[2]
  if (length(aux) == 1) {
    f <- NULL
  } else {
    f <- data[, namef]
  }
  newdata <- data
  data <- na.omit(data[, c(ffr$response, varnames)])
  newdata <- na.omit(newdata[, varnames])
  n <- nrow(data)
  
  #if (is.null(seed)) {
  #  set.seed(NULL)
  #  seed <- .Random.seed[3]
  #}
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  umatrix <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
  
  if (is.null(f)) 
    f <- rep(1, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  
  res <- list()
  
  for (i in etiquetas) {
    yy <- data[, 1][f == i]
    xx <- data[, 2][f == i]
    n <- length(xx)
    w <- rep(1, n)
    fit <- .Fortran("allotest_", 
                    x = as.double(xx), 
                    y = as.double(yy), 
                    w = as.double(w), 
                    n = as.integer(n), 
                    kbin = as.integer(kbin), 
                    nboot = as.integer(nboot), 
                    seed = as.integer(seed), 
                    T = as.double(-1), 
                    pvalue = as.double(-1),
                    umatrix = as.double(umatrix)
                    )
    
    res[[i]] <- list(statistic = fit$T, pvalue = fit$pvalue)
  }
  
  
  est <- c()
  p <- c()
  for (i in 1:length(res)) {
    est[i] <- res[[i]]$statistic
    p[i] <- res[[i]]$pvalue
  }
  
  kk <- cbind(est, p)
  result <- matrix(kk, ncol = 2, nrow = length(res))
  colnames(result) <- c("Statistic", "pvalue")
  
  
  # factores=paste('Factor',1:length(res))
  factores <- paste("Level", etiquetas[1:length(etiquetas)])
  
  
  for (i in 1:length(res)) {
    rownames(result) <- c(factores)
  }
  
  return(result)
  
} 