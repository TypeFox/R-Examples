##' Compare two link functions for the binomial distribution.
##'
##' Finds the closest (after suitable scaling) of the link function
##' \code{d2} to \code{d1}. If \code{prob1} is provided, then this
##' corresponds to probabilities under link \code{d1} and the
##' probabilities are transformed to \code{d2}. Also creates a plot of
##' the two link functions.
##' @title Compare two link functions for the binomial distribution
##' @param d1 Either a positive number of one of "logit" or "probit"
##' defining the link function. If positive number, this corresponds
##' to the \code{robit(d1)} link.
##' @param d2 Same as \code{d1}.
##' @param a Beginning of range of points to evaluate (and plot) the
##' two link functions.
##' @param b End of range of points to evaluate (and plot) the two
##' link functions.
##' @param n Number of points for evaluating (and plot) the two link
##' functions.
##' @param prob1 Binomial probabilities corresponding to the first
##' link function.
##' @param plot Whether the two link functions should be plotted. If
##' so, two plots are created: the left plot shows two curves
##' corresponding to the two link functions (black is the first) and
##' the right plot shows their difference (first - second).
##' @return A list with the following elements
##' \itemize{
##'   \item \code{scale} The optimal scaling \code{c} such that
##'     \code{max(abs(link(x/c,d1) - link(x,d2)))} is minimised for
##'     \code{x = seq(a,b,length.out=n)}
##'   \item \code{maxdiff} The maximum difference between the two links.
##'   \item \code{prob2} The corresponding probabilities from \code{prob1}
##'     to the second link.}
##' @examples
##' \dontrun{
##' comparebinlinks("logit", 7) # The robit(7) approximates logit
##' comparebinlinks("probit", 1, prob1 = c(.5, .6, .7, .8, .9))
##' }
##' @importFrom stats plogis pnorm pt qlogis qnorm qt optimize
##' @importFrom graphics curve par
##' @export 
comparebinlinks <- function (d1,d2,a=-8,b=8,n=2001,prob1=NULL,plot=TRUE) {
  ## d1, d2 - Degrees of freedom
  ## a, b - Interval
  ## n - Number of points
  link <- function(x,s,d) {
    if (d < 0) {
      plogis(x,,s)
    } else if (d == 0) {
      pnorm(x,,s)
    } else {
      pt(x/s,d)
    }
  }
  linkinv <- function(p,s,d) {
    if (d < 0) {
      qlogis(p,0,s)
    } else if (d == 0) {
      qnorm(p,,s)
    } else {
      qt(p,d)*s
    }
  }
  link2num <- function (d) {
    if (length(d) != 1) stop('Link input must be a numeric or character
of lenght 1')
    links <- c('logit','probit')
    if (is.character(d) | is.factor(d)) {
      dnm <- match.arg(d,links)
      d <- which(dnm == links) - length(links)
    } else if (d <= 0) {
      stop (paste('Link function must be either a positive number or',
                  paste(paste('"',links,'"',sep=''),collapse=', ')))
    } else {
      dnm <- paste('robit(',d,')',sep='')
    }
    attr(d,'name') <- dnm
    d
  }
  d1 <- link2num(d1)
  d2 <- link2num(d2)
  if (d1 == d2) {
    out <- list()
    out$scale <- 1
    out$maxdiff <- 0
    out$prob2 <- prob1
    if (plot) {
      oldpar <- par(mfrow=c(1,2))
      on.exit(par(oldpar), add = TRUE)
      curve(link(x,1,d1),a,b,n,col=2,xlab='Linear predictor',ylab='',main='Link')
      curve(0*x,a,b,n=2,col=2,xlab='Linear predictor',ylab='',main='Difference')
    }
  } else {
    x <- seq(a,b,length=n)
    f <- function(s) max(abs(link(x,1,d2)-link(x,s,d1)))
    fop <- optimize(f,c(1e-4,1e4))
    if (plot) {
      oldpar <- par(mfrow=c(1,2))
      on.exit(par(oldpar), add = TRUE)
      curve(link(x,fop$minimum,d1),a,b,n,
            xlab='Linear predictor',ylab='',main='Links')
      curve(link(x,1,d2),n=n,add=TRUE,col=2)
      curve(link(x,1,d2)-link(x,fop$minimum,d1),a,b,n,col=2,
            xlab='Linear predictor',ylab='',main='Difference')
      out <- list()
      out$scale <- fop$minimum
      out$maxdiff <- fop$objective
      if (is.null(prob1)) {
        out$prob2 <- NULL
      } else {
        out$prob2 <- link(linkinv(prob1,fop$minimum,d1),1,d2)
      }
    }
  }
  out
}
