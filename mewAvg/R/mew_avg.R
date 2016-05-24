## This function packages the MEW process into a single function as
## opposed to requiring the user to call the function mewInit,
## mewAccum and mewMean themselves, appropriately.  The downfall of
## this interface is that the user cannot run the algorithm for some
## number of iterations, pause, assess convergence of the mean and
## then pick up from where they paused.  To accomplish that see the
## examples associated with the mewMean function.

#' @title Convenience wrapper for the MEW process
#'
#' @description Packages the process of calling \code{mewInit},
#' looping through the random vectors calling \code{mewAccum} for each
#' one and calling \code{mewMean} when desired.
#'
#' @details The function \code{f} should generate the sequence of
#' random vectors one at a time.  The returned value from a single call
#' should be a list with at least one element.  The first element
#' should be a numeric vector of length \code{n.xx} (the next vector
#' in the sequence), and the remaining elements should be the updated
#' arguments for the next call to \code{f}, named appropriately for
#' the argument of \code{f} to update.  The 'Examples' section
#' provides further guidance.
#'
#' The downfall of this interface is that the user cannot run the
#' algorithm for some number of iterations, pause, assess convergence
#' of the mean and then pick up from where they paused.  To accomplish
#' that see the examples associated with the \code{mewMean} function.
#'
#' @param f (function) A user defined R function.  See the 'Details'
#' section for more on defining this function
#'
#' @param n.bin (scalar integer) The fixed number of bins to use to
#' define the moving expanding window
#'
#' @param n.xx (scalar integer) The length of the numeric vector
#' returned by \code{f}
#'
#' @param ff (scalar double) The fraction of the samples to included
#' in each window
#'
#' @param n.save (scalar integer OR NULL)The number of estimates to
#' save and return.  The default value is NULL since this argument can
#' be derived from \code{i.to.save}.  The argument is kept for
#' compatibility with older versions of this package
#'
#' @param n.iter (scalar integer OR NULL) The number of times to call
#' \code{f}.  The default value is NULL since this argument can be
#' derived from \code{i.to.save}.  The argument is kept for
#' compatibility with older versions of this package
#'
#' @param i.to.save (vector integer length n.iter) A vector of zeros
#' and ones of length \code{n.iter} where position \code{i} is 1 if an
#' average should be calculated and saved at iteration i, and zero
#' otherwise
#'
#' @param ... The initial named arguments to \code{f}.
#'
#' @return A matrix of dimension \code{n.save} by \code{n.xx}
#' containing the saved averages
#'
#' @examples
#' MyFun <- function (k) {
#'
#'  value <- runif(n=2)
#'  value[1] <- ((cos(value[1]*2*pi))^2)*(1 - exp(-0.01*k))
#'  value[2] <- (-((sin(value[2]*2*pi))^2))*(1 - exp(-0.01*k))
#'
#'  k <- k + 1
#'
#'  return(list(value=value, k=k))
#' }
#'
#' i.to.save <- seq(from=1, to=1025, by=32)
#' tmp <- rep(x=0, times=1025)
#' tmp[i.to.save] <- 1
#' i.to.save <- tmp
#'
#' mean.vals <- mewAvg(f=MyFun,
#'                     n.bin=4,
#'                     n.xx=2,
#'                     ff=0.5,
#'                     n.save=sum(i.to.save),
#'                     n.iter=length(i.to.save),
#'                     i.to.save=i.to.save,
#'                     k=1)
#'
#' plot(c(1:sum(i.to.save),
#'        1:sum(i.to.save)),
#'      c(mean.vals[, 1],
#'        mean.vals[, 2]),
#'      type="n",
#'      xlab="Saved Iter",
#'      ylab="Mean")
#' points(1:sum(i.to.save),
#'        mean.vals[, 1])
#' points(1:sum(i.to.save),
#'        mean.vals[, 2])
#'
#' ## an AR(1) process
#'
#' ArOne <- function (x.old, phi, sig.eps) {
#'
#'   value <- phi*x.old + rnorm(n=1, mean=0, sd=sig.eps)
#'
#'   return(list(value=value, x.old=value))
#' }
#'
#' mean.vals.ar1 <- mewAvg(f=ArOne,
#'                         n.bin=4,
#'                         n.xx=1,
#'                         ff=0.5,
#'                         n.save=sum(i.to.save),
#'                         n.iter=length(i.to.save),
#'                         i.to.save=i.to.save,
#'                         x.old=0,
#'                         phi=0.5,
#'                         sig.eps=1)
#'
#' plot(x=c(1, sum(i.to.save)),
#'      y=c(-0.5, 0.5),
#'      xlab="Saved Iter",
#'      ylab="Mean",
#'      type="n")
#' points(x=1:sum(i.to.save),
#'        y=mean.vals.ar1)
#' abline(h=0, col="red")
#'
#' @export
mewAvg <- function (f, n.bin, n.xx, ff,
                    n.save = NULL, n.iter = NULL,
                    i.to.save, ...) {

  if (is.null(n.save)) {

    n_save <- sum(i.to.save)
  } else {

    n_save <- n.save
  }

  if (is.null(n.iter)) {

    n_iter <- length(i.to.save)
  } else {

    n_iter <- n.iter
  }

  el_args <- list(...)

  if (length(el_args) != length(formals(f))) {

    stop(paste0("There must be the same number of arguments\n",
                "supplied by ... as there are for f"))
  }

  av <- mewInit(n_bin = n.bin,
                n_xx = n.xx,
                ff = ff)

  results <- matrix(double(n_save*n.xx),
                    nrow = n_save,
                    ncol = n.xx)

  rho <- new.env()

  n_complete <- 0

  n_already_saved <- 0

  while (n_complete < n_iter) {

    n_complete <- n_complete + 1

    if (length(el_args) > 0){

      for (i in 1:length(el_args)) {

        assign(x = names(el_args)[i],
               value = el_args[[i]],
               envir = rho)
      }
    }

    f_eval <- eval(parse(text = paste0("f(",
                             paste(names(formals(f)), collapse = ","),
                             ")")),
                   envir = rho)

    xx <- f_eval[[1]]

    av <- mewAccum(xx = xx, av = av)

    if (i.to.save[n_complete] == 1) {

      n_already_saved <- n_already_saved + 1

      av <- mewMean(av)

      results[n_already_saved, ] <- mewGetMean(av)
    }

    el_args <- NULL
    if (length(f_eval) > 1) {

      for (i in 2:length(f_eval)) {

        el_args <- c(el_args, f_eval[i])
      }
    }
  }

  return(results)
}
