#' @name LambertW-toolkit
#' @title Do-it-yourself toolkit for Lambert W\eqn{ \times} F distribution
#' @description
#'
#' \strong{IMPORTANT:} This toolkit functionality is still under active
#'     development; function names, arguments, return values, etc. may change.
#' 
#' This do-it-yourself Lambert W\eqn{ \times} F toolkit implements the flexible
#'     input/output framework of Lambert W \eqn{\times} F random variables (see
#'     References).  Using a modular approach, it allows users to create their
#'     own Lambert W \eqn{\times} 'MyFavoriteDistribution' RVs. See Details
#'     below.
#' 
#' If the distribution you inted to use is not already implemented
#' (\code{\link{get_distnames}}), then you can create it:
#' \describe{
#' \item{create input:}{use \code{\link{create_LambertW_input}} with your
#' favorite distribution,}
#' \item{create output:}{pass it as an input argument to \code{\link{create_LambertW_output}},}
#' \item{use output:}{use Rs standard functionality for distributions
#' such as random number generation (\code{rY}), pdf (\code{dY}) and cdf
#' (\code{pY}), quantile function (\code{qY}), etc. for this newly generated
#' Lambert W \eqn{\times} 'MyFavoriteDistribution'.}
#' }
#' 
#' @details
#' \code{\link{create_LambertW_output}} takes an object of class
#'     \code{LambertW_input} and creates a class \code{LambertW_output} for
#'     standard distributions as well as the user-defined distribution. This
#'     \code{LambertW_output} represents the RV Y \eqn{\sim} Lambert W
#'     \eqn{\times} 'MyFavoriteDistribution' with all its properties and R
#'     functionality, such as random number generation (\code{rY}), pdf
#'     (\code{dY}) and cdf (\code{pY}), etc.
#' 
#' @inheritParams common-arguments
#' @author Georg M. Goerg
#' @keywords univar distribution datagen models
#' @examples
#' 
#' # create a Gaussian N(1, 2) input
#' Gauss.input <- create_LambertW_input("normal", beta = c(1, 2))
#' 
#' # create a heavy-tailed version of a normal
#' # gamma = 0, alpha = 1 are set by default; beta comes from input
#' params <- list(delta = c(0.3)) 
#' LW.Gauss <- create_LambertW_output(LambertW.input = Gauss.input, 
#'                                    theta = params)
#' LW.Gauss
#' 
#' op <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 1), mar = c(3, 3, 2, 1))
#' curve(LW.Gauss$d(x, params), -7, 10, col = "red")
#' # parameter will get detected automatically from the input
#' curve(LW.Gauss$d(x), -7, 10, col = "blue") # same in blue; 
#' 
#' # compare to the input case (i.e. set delta = 0)
#' params.0 <- params 
#' params.0$delta <- 0
#' 
#' # to evaluate the RV at a different parameter value, 
#' # it is necessary to pass the new parameter
#' curve(LW.Gauss$d(x, params.0), -7, 10, add = TRUE, col = 1) #' par(op)
#' 
#' curve(LW.Gauss$p(x, params), -7, 10, col = "red")
#' curve(LW.Gauss$p(x, params.0), -7, 10, add = TRUE, col = 1)
#' 
#' test_normality(LW.Gauss$r(n = 100), add.legend = FALSE)
#' 
#' ## generate a positively skewed version of a shifted, scaled t_3
#' t.input <- create_LambertW_input("t", beta = c(2, 1, 3))
#' t.input
#' params <- list(gamma = 0.05) # skew it
#' LW.t <- create_LambertW_output(LambertW.input = t.input, theta = params)
#' LW.t
#' 
#' plot(t.input$d, -7, 11, col = 1)
#' plot(LW.t$d, -7, 11, col = 2, add = TRUE)
#' abline(v = t.input$beta["location"], lty = 2)
#' 
#' # draw samples from the skewed t_3
#' yy <- LW.t$r(n = 100)
#' test_normality(yy)
#' 
#' ### create a skewed exponential distribution
#' exp.input <- create_LambertW_input("exp", beta = 1)
#' plot(exp.input)
#' params <- list(gamma = 0.2)
#' LW.exp <- create_LambertW_output(exp.input, theta = params)
#' plot(LW.exp)
#' 
#' # create a heavy-tail exponential distribution
#' params <- list(delta = 0.2)
#' LW.exp <- create_LambertW_output(exp.input, theta = params)
#' plot(LW.exp)
#' 
#' # create a skewed chi-square distribution with 5 df
#' chi.input <- create_LambertW_input("chisq", beta = 5)
#' plot(chi.input)
#' params <- list(gamma = sqrt(2)*0.2)
#' LW.chi <- create_LambertW_output(chi.input, theta = params)
#' plot(LW.chi)
#' 
#' 
NULL
