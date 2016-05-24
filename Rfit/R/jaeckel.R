#' Function to Minimize Jaeckel's Dispersion Function
#' 
#' Uses the built-in function \code{optim} to minimize Jaeckel's dispersion
#' function.  Alternates between \code{CG} and \code{BFGS} steps and initially
#' takes a number of IRLS steps.
#' 
#' This function is meant to mimic the minization algorithm implemented in RGLM
#' (Kapenga, et. al. 1988). The main loop of the function alternates between CG
#' and BFGS estimation methods.  To speed up convergence, it first takes a
#' number of iterated reweighted least squares (IRLS) steps which generally
#' gets close to the solution in a small number of steps.  Using IRLS for rank
#' regression was first considered by Cheng and Hettmansperger (1983).  See
#' also Sievers and Abebe (2004).
#' 
#' @param x n by p design matrix
#' @param y n by 1 response vector
#' @param beta0 intial estimate
#' @param scores object of class 'scores'
#' @return Results of \code{optim} are returned.
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{optim}}, \code{\link{rfit}}
#' @references Cheng, K. S. and Hettmansperger, T. P. (1983), Weighted
#' Least-Squares Rank Regression, \emph{Communications in Statistics, Part A -
#' Theory and Methods}, 12, 1069-1086.
#' 
#' Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric
#' Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Jaeckel, L. A. (1972), Estimating regression coefficients by minimizing the
#' dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449
#' - 1458.
#' 
#' Kapenga, J. A., McKean, J. W., and Vidmar, T. J. (1988), \emph{RGLM: Users
#' Manual}, Statist. Assoc. Short Course on Robust Statistical Procedures for
#' the Analysis of Linear and Nonlinear Models, New Orleans.
#' 
#' Sievers, J. and Abebe, A. (2004), Rank Estimation of Regression Coefficients
#' Using Iterated Reweighted Least Squares, \emph{Journal of Statistical
#' Computation and Simulation}, 74, 821-831.
#' @examples
#' 
#' ##  This is a internal function.  See rfit for user-level examples.
#' 
#' ## The function is currently defined as
#' function (x, y, beta0 = rq(y ~ x - 1)$coef, scores = Rfit::wscores, 
#'     maxiter = 100, irls0 = 10, BFGS0 = 20, stepCG = 5, stepBFGS = 2) 
#' {
#'     x <- x - outer(rep(1, nrow(x)), apply(x, 2, mean))
#'     beta0 <- irls(x, y, beta0, max.iter = irls0)
#'     if (BFGS0 < 1) 
#'         BFGS0 <- 1
#'     fit <- optim(beta0, disp, method = "BFGS", x = x, y = y, 
#'         scores = scores, control = list(maxit = BFGS0))
#'     iter <- 0
#'     while (fit$convergence && iter < maxiter) {
#'         iter <- iter + 1
#'         fit <- optim(fit$par, disp, method = "CG", x = x, y = y, 
#'             scores = scores, control = list(maxit = stepCG))
#'         fit <- optim(fit$par, disp, method = "BFGS", x = x, y = y, 
#'             scores = scores, control = list(maxit = stepBFGS))
#'     }
#'     optim(fit$par, disp, method = "BFGS", x = x, y = y, scores = scores)
#'   }
#' 
#' @export jaeckel
jaeckel <- function (x, y, beta0 = rq(y ~ x)$coef[2:(ncol(x) + 1)], 
  scores = Rfit::wscores, ...) {

  scrs <- getScores(scores,seq_len(length(y))/(length(y)+1))

  j.grad <- function (x, y, beta,scores,...) {
    x <- as.matrix(x)
    e <- y - x %*% beta
    r <- rank(e, ties.method = "first")/(length(e) + 1)
    crossprod(x,-1*getScores(scores,r) )
  }


  j.disp <- function (beta, x, y,scrs,...) {
    e <- y - x %*% beta
    drop(crossprod(e[order(e)],scrs))
  }

  sd.y <- sd(y)
  ystar <- y/sd.y

  fit0 <- optim(beta0/sd.y, j.disp, method = "BFGS", x = x, y = ystar, scrs=scrs,
    gr=j.grad,scores=scores,...)

  optim(fit0$par*sd.y, j.disp, method = "BFGS", x = x, y = y, scrs=scrs,
    gr=j.grad,scores=scores,...)

}

