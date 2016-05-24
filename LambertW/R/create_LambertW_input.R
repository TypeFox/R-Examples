#' @rdname LambertW-toolkit
#' 
#' @details
#' 
#' \code{create_LambertW_input} allows users to define their own Lambert
#'     W\eqn{\times} F distribution by supplying the necessary functions about
#'     the input random variable \eqn{U} and \eqn{\boldsymbol \beta}.  Here
#'     \eqn{U} is the zero mean and/or unit variance version of \eqn{X \sim
#'     F_X(x \mid \boldsymbol \beta)} (see References).
#' 
#' The argument \code{input.u} must be a list containing all of the following:
#' \describe{
#' \item{\code{beta2tau}}{ R function of \code{(beta)}: converts \eqn{\boldsymbol \beta} to \eqn{\tau} for the
#' user defined distribution }
#' \item{\code{distname}}{ optional; users can specify the name 
#' of their input distribution. By default it's called \code{"MyFavoriteDistribution"}. 
#' The distribution name will be used in plots and summaries of the Lambert W\eqn{\times} F 
#' input (and output) object.}
#' \item{\code{is.non.negative}}{ logical; users should specify whether the
#' distribution is for non-negative random variables or not.  This will help
#' for plotting and theoretical quantile computation.}
#' \item{\code{d}}{ R function of \code{(u, beta)}: probability density function (pdf) of U,}
#' \item{\code{p}}{ R function of \code{(u, beta)}: cumulative distribution function (cdf) of U,}
#' \item{\code{q}}{ R function of \code{(p, beta)}: quantile function of U,}
#' \item{\code{r}}{ R function \code{(n, beta)}: random number generator for U,}
#' }
#' 
#' @inheritParams common-arguments
#' 
#' @param input.u optional; users can make their own 'Lambert W x F'
#'     distribution by supplying the necessary functions. See Description for
#'     details.

#' @keywords univar distribution datagen models
#' @export
#' @examples
#' # a demo on how a user-defined U input needs to look like
#' user.tmp <- list(d = function(u, beta) dnorm(u),
#'                  r = function(n, beta) rnorm(n),
#'                  p = function(u, beta) pnorm(u),
#'                  q = function(p, beta) qnorm(p),
#'                  beta2tau = function(beta) {
#'                    c(mu_x = beta[1], sigma_x = beta[2], 
#'                      gamma = 0, alpha = 1, delta = 0)
#'                    },
#'                  distname = "MyNormal",
#'                  is.non.negative = FALSE)
#' my.input <- create_LambertW_input(input.u = user.tmp, beta = c(0, 1))
#' my.input
#' plot(my.input)
#' 

create_LambertW_input <- function(distname = NULL, beta, 
                                  input.u = list(beta2tau = NULL, 
                                                 d = NULL, p = NULL, r = NULL, q = NULL,
                                                 distname = "MyFavoriteDistribution",
                                                 is.non.negative = FALSE)) {
  
  required.input.names <- c("beta2tau",
                            "d", "q", "r", "q",
                            "distname",
                            "is.non.negative")
  stopifnot(all(match(required.input.names, names(input.u))))
  if (is.null(input.u$distname)) {
    input.u[["distname"]] <- "MyFavoriteDistribution"
  }
  all.U.fcts.available <- all(!sapply(input.u, is.null))
  
  if (!all.U.fcts.available && is.null(distname)) {
    stop("You must either provide a distribution by name or provide your own user defined functions",
         " to that define random variable U and the beta2tau transformation.")    
  } else if (!is.null(distname) && all.U.fcts.available) {
    stop("Please choose a distribution _either_ by name or provide your own functions for U.", 
         " Do not provie both.")
  }
  
  user.defined <- FALSE
  if (all.U.fcts.available) {
    user.defined <- TRUE
    # specify local U fcts
    dU_tmp <- function(u) input.u$d(u, beta = beta)
    rU_tmp <- function(u) input.u$r(u, beta = beta)
    pU_tmp <- function(u) input.u$p(u, beta = beta)
    qU_tmp <- function(u) input.u$q(u, beta = beta)
    beta2tau_tmp <- function(beta) input.u$beta2tau(beta = beta)
    is.non.negative <- input.u$is.non.negative
    if (any(names(beta))) {
      warning("Your 'beta' vector does not have proper 'names'. Please fix.")
    }
  } else if (!is.null(distname)) {
    dist.family <- get_distname_family(distname)
    is.non.negative <- dist.family$is.non.negative
    dU_tmp <- function(u) dU(u, beta = beta, distname = distname)
    rU_tmp <- function(u) rU(u, beta = beta, distname = distname)
    pU_tmp <- function(u) pU(u, beta = beta, distname = distname)
    qU_tmp <- function(u) qU(u, beta = beta, distname = distname)
    beta2tau_tmp <- function(beta) beta2tau(beta = beta, distname = distname)
    names(beta) <- get_beta_names(distname)
    
  } else {
    stop("Something went wrong in create_LambertW_input(). Either 'distname' or user-defined functions",
         " are incorrectly specified.")
  }

  result <- list(beta = beta,
                 tau = beta2tau_tmp(beta = beta),
                 distname = ifelse(is.null(distname), input.u[["distname"]], distname),
                 user.defined = user.defined,
                 is.non.negative = is.non.negative)

  rX <- function(n, beta = result$beta) {
    tau <- beta2tau_tmp(beta = beta)
    sim.vals <- rU_tmp(n) * tau["sigma_x"] + tau["mu_x"]
    names(sim.vals) <- NULL
    return(sim.vals)
  }
  pX <- function(x, beta = result$beta) {
    tau <- beta2tau_tmp(beta = beta)
    cdf.val <- pU_tmp((x - tau["mu_x"]) / tau["sigma_x"])
    names(cdf.val) <- NULL
    return(cdf.val)
  }
  dX <- function(x, beta = result$beta) {
    tau <- beta2tau_tmp(beta = beta)
    pdf.val <- dU_tmp((x - tau["mu_x"]) / tau["sigma_x"])
    names(pdf.val) <- NULL
    return(pdf.val)
  }
  qX <- function(p, beta = result$beta) {
    tau <- beta2tau_tmp(beta = beta)
    quant.val <- qU_tmp(p) * tau["sigma_x"] + tau["mu_x"]
    names(quant.val) <- NULL
    return(quant.val)
  }
  
  # add to output
  result <- c(result,
              list(U = list(d = dU_tmp,
                            p = pU_tmp,
                            r = rU_tmp,
                            q = qU_tmp)))
  result <- c(result,
              list(d = dX,
                   p = pX,
                   r = rX,
                   q = qX))

  result <- c(result,
              list(distname.with.beta = paste0(result$distname, "(", 
                                               paste(round(result$beta, 2), collapse = ","), ")"),
                   beta2tau = beta2tau_tmp))
  class(result) <- "LambertW_input"
  invisible(result)
}
