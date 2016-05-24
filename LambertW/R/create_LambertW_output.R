#' @rdname LambertW-toolkit
#' @description
#' 
#' \code{create_LambertW_output} converts the input \code{LambertW_input}
#'     representing random variable \eqn{X \sim F_X} to the Lambert W
#'     \eqn{\times} \eqn{F_X} output.
#' 
#' @inheritParams common-arguments
#' @param LambertW.input an object of class \code{LambertW_input}
#' @return 
#' \code{create_LambertW_output} returns a list of class \code{LambertW_output} 
#' with values that are (for the most part) functions themselves (see Examples): 
#' 
#' \item{d}{ pdf of Y \eqn{\sim} Lambert W \eqn{\times} 'MyFavoriteDistribution',} 
#' \item{p}{ cdf of Y,} 
#' \item{q}{ quantile function for Y,} 
#' \item{r}{ random number generator for Y,}
#' \item{distname}{ character string with the name of the new distribution.
#' Format: "Lambert W x 'MyFavoriteDistribution'",} 
#' \item{beta, theta}{see Arguments,}
#' \item{distname.with.beta}{name of the new distribution
#' including the parameter \code{beta}. Format: "Lambert W x 'MyFavoriteDistribution'(beta)".}
#' 
#' @keywords univar distribution datagen models
#' @export
#' 

create_LambertW_output <- function(LambertW.input = NULL, theta = NULL, 
                                   distname = LambertW.input$distname) {
  
  if (is.null(LambertW.input) && is.null(distname)) {
    stop("'LambertW.input' is missing. \n Please create an object of class ",
         " 'LambertW_input' first (using 'create_LambertW_input()` and then ",
         "pass it as \n ",
         "'create_LambertW_output(LambertW.input = my_LambertW_input, theta = ...).")
  } 
  
  if (!LambertW.input$user.defined) {
    check_distname(distname)
  }
  
  if (is.null(theta)) {
    cat("You have not specified the theta of the Lambert W RV transformation.",
        "\n By default your LambertW output will be identical to the ",
        "LambertW input.")
    theta <- list(beta = LambertW.input$beta)
  }
  theta <- complete_theta(theta, LambertW.input = LambertW.input)

  if (!is.null(LambertW.input)) {
    distname <- LambertW.input$distname
    check_theta(theta = theta, distname = distname)
  }
  
  if (is.null(LambertW.input)) {
    LambertW.input <- create_LambertW_input(distname = distname,
                                            beta = theta$beta)
  }
  
  tau.tmp <- LambertW.input$beta2tau(LambertW.input$beta)
  # update with theta parameters
  theta.tmp <- theta
  # pretend its a normal distribution, and then convert to tau again
  theta.tmp$beta <- tau.tmp[c("mu_x", "sigma_x")]
  tau.update <- theta2tau(theta.tmp, distname = "normal")

  results <- list(call = match.call(),
                  theta = theta,
                  beta = LambertW.input$beta,
                  tau = tau.update,
                  input.distname = LambertW.input$distname,
                  distname = paste("Lambert W x", LambertW.input$distname))
  results$type <- tau2type(results$tau)
  results$distname.with.beta <- paste0(results$distname, "(", 
                                       paste(round(results$beta, 2),
                                             collapse = ","), ")")
  
  rY <- function(n, theta = results$theta) {
    theta <- complete_theta(theta, LambertW.input = LambertW.input)
    tau.tmp <- LambertW.input$beta2tau(theta$beta)
    # update with theta parameters
    theta.tmp <- theta
    # pretend its a normal distribution, and then convert to tau again
    theta.tmp$beta <- tau.tmp[c("mu_x", "sigma_x")]
    tau.update <- theta2tau(theta.tmp, distname = "normal")
    sim.vals <- rLambertW(n = n, theta = theta, input.u = LambertW.input$U$r(n),
                          tau = tau.update)
    names(sim.vals) <- NULL
    return(sim.vals)
  }
  
  dY <- function(y, theta = results$theta) {
    theta <- complete_theta(theta, LambertW.input = LambertW.input)
    
    tau.tmp <- LambertW.input$beta2tau(theta$beta)
    # update with theta parameters
    theta.tmp <- theta
    # pretend its a normal distribution, and then convert to tau again
    theta.tmp$beta <- tau.tmp[c("mu_x", "sigma_x")]
    tau.update <- theta2tau(theta.tmp, distname = "normal")
    pdf.vals <- dLambertW(y, theta = theta, input.u = LambertW.input$U$d,
                          tau = tau.update)
    names(pdf.vals) <- NULL
    return(pdf.vals)
  }
  
  pY <- function(y, theta = results$theta) {
    theta <- complete_theta(theta, LambertW.input = LambertW.input)
    tau.tmp <- LambertW.input$beta2tau(theta$beta)
    # update with theta parameters
    theta.tmp <- theta
    # pretend its a normal distribution, and then convert to tau again
    theta.tmp$beta <- tau.tmp[c("mu_x", "sigma_x")]
    tau.update <- theta2tau(theta.tmp, distname = "normal")
    cdf.vals <- pLambertW(y, theta = theta, input.u = LambertW.input$U$p,
                          tau = tau.update)
    names(cdf.vals) <- NULL
    return(cdf.vals)
  }
  
  qY <- function(p, theta = results$theta) {
    theta <- complete_theta(theta, LambertW.input = LambertW.input)
    
    tau.tmp <- LambertW.input$beta2tau(theta$beta)
    # update with theta parameters
    theta.tmp <- theta
    # pretend its a normal distribution, and then convert to tau again
    theta.tmp$beta <- tau.tmp[c("mu_x", "sigma_x")]
    tau.update <- theta2tau(theta.tmp, distname = "normal")
    quant.vals <- qLambertW(p, theta = theta, input.u = LambertW.input$U$q,
                            tau = tau.update, 
                            is.non.negative = LambertW.input$is.non.negative)
    names(quant.vals)
    return(quant.vals)
  }
  
  results <- c(results, 
               list(d = dY,
                    p = pY,
                    r = rY,
                    q = qY))  
  class(results) <- "LambertW_output"
  invisible(results)
}
