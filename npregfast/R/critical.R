#' Critical points of the regression function
#' 
#' @description This function draws inference about some critical point in 
#' the support of \eqn{X} which is associated  with  some features of the regression
#' function (e.g., minimum, maximum or inflection points which indicate changes
#' in the sign of curvature). Returns the value of the covariate \code{x} 
#' which maximizes the estimate of the function, the value of the covariate 
#' \code{x} which maximizes the first derivative and the value of the covariate 
#' \code{x} which equals the second derivative to zero, for each level of the 
#' factor.
#'
#'@param model Parametric or nonparametric regression out 
#' obtained by \code{\link{frfast}} function.
#'@param der Number which determines any inference process. By default
#' \code{der} is \code{NULL}. If this term is \code{0}, the calculation is for
#' the point which maximize the estimate. If it is \code{1} it is 
#' designed for the first derivative and if it is \code{2}, it returns the point
#'  which equals the second derivative to zero.
#' 
#'@return An object is returned with the following elements:
#'\item{Estimation}{ \code{x} value which maximize the regression function with 
#'their 95\% confidence intervals (for each level).}
#'\item{First_der}{\code{x} value which maximize the first derivative with their
#' 95\% confidence intervals (for each level).}
#' \item{Second_der}{\code{x} value which equals the second derivative to zero 
#' with their 95\% confidence intervals (for each level).}
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'
#' @references 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#'
#'@examples
#'library(npregfast)
#'data(barnacle)
#'
#'fit <- frfast(DW ~ RC, data = barnacle) # without interactions
#'critical(fit)
#'critical(fit, der = 0)
#'critical(fit, der = 1)
#'critical(fit, der = 2)
#'
#' # fit2 <- frfast(DW ~ RC : F, data = barnacle) # with interactions
#' # critical(fit2)
#' # critical(fit2, der = 0)
#' # critical(fit2, der = 1)
#' # critical(fit2, der = 2)
#'
#'@export




critical <- function(model, der = NULL) {
  
  if(length(der) > 1){
    stop("Argument \"der\" have to be a length-one vector")
  }
  
  if(!is.null(der) & !isTRUE(der %in% c(0, 1, 2))) {
    stop("",paste(der)," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  }
  
  nf <- model$nf
  jnf <- c()
  model$max[model$max == 9999] <- NA
  model$maxl[model$maxl == 9999] <- NA
  model$maxu[model$maxu == 9999] <- NA
  factores <- c()
  jnf <- c()
  if (is.null(der)) {
    res <- matrix(ncol = 3, nrow = nf)
    k <- 1
    for (i in 1:nf) {
      res[i, 1] <- c(model$max[k, i])
      res[i, 2] <- c(model$maxl[k, i])
      res[i, 3] <- c(model$maxu[k, i])
    }
    res2 <- matrix(ncol = 3, nrow = nf)
    k <- 2
    for (i in 1:nf) {
      res2[i, 1] <- c(model$max[k, i])
      res2[i, 2] <- c(model$maxl[k, i])
      res2[i, 3] <- c(model$maxu[k, i])
    }
    res3 <- matrix(ncol = 3, nrow = nf)
    k <- 3
    
    for (i in 1:nf) {
      res3[i, 1] <- c(model$max[k, i])
      res3[i, 2] <- c(model$maxl[k, i])
      res3[i, 3] <- c(model$maxu[k, i])
      jnf[i] <- which(model$label == model$label[i])
      factores[i] <- paste("Level", model$label[jnf[i]])
    }
    
    colnames(res) <- c("Critical", "Lwr", "Upr")
    colnames(res2) <- c("Critical", "Lwr", "Upr")
    colnames(res3) <- c("Critical", "Lwr", "Upr")
    rownames(res) <- c(factores)
    rownames(res2) <- c(factores)
    rownames(res3) <- c(factores)
    return(list(Estimation = res, First_der = res2, Second_der = res3))
  } else {
    der <- der + 1
    res <- matrix(ncol = 3, nrow = nf * length(der))
    k <- der
    a <- 1
    
    for (j in der) {
      for (i in 1:nf) {
        if (a == 2) {
          ii <- nf + i
        } else {
          ii <- i
        }
        res[ii, 1] <- c(model$max[j, i])
        res[ii, 2] <- c(model$maxl[j, i])
        res[ii, 3] <- c(model$maxu[j, i])
        jnf[i] <- which(model$label == model$label[i])
        factores[i] <- paste("Level", model$label[jnf[i]])
      }
      a <- 2
    }
    colnames(res) <- c("Critical", "Lwr", "Upr")
    rownames(res) <- c(rep(factores, length(der)))
    return(res)
  }
}
