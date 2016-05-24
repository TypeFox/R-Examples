#' @title Calculate option Greeks
#'
#' @description The functions \code{greeks} and \code{greeks2} provide
#' two different calling conventions for computing a full set of
#' option Greeks. \code{greeks} simply requires entering a pricing function
#' with parameters. \code{greeks2} requires the use of named parameter
#' entries. The function \code{bsopt} calls \code{greeks2} to
#' produce a full set of prices and greeks for calls and puts. These
#' functions are all vectorized, the only restriction being that the
#' functions will produce an error if the recycling rule can not be
#' used safely (that is, if parameter vector lengths are not integer
#' multiples of one another).
#'
#' @name greeks
#' @aliases bsopt greeks greeks2
#'
#' @return A named list of Black-Scholes option prices and Greeks.
#'
#' @details Numerical derivatives are calculated using a simple
#' difference. This can create numerical problems in edge cases. It
#' might be good to use the package numDeriv or some other more
#' sophisticated calculation, but the current approach works well with
#' vectorization.
#' 
#' @usage
#' bsopt(s, k, v, r, tt, d)
#' greeks(f)
#' # must used named list entries:
#' greeks2(fn, ...)
#'
#' @param s Stock price
#' @param k Strike price of the option
#' @param v Volatility of the stock, defined as the annualized
#'     standard deviation of the continuously-compounded return
#' @param r Annual continuously-compounded risk-free interest rate
#' @param tt Time to maturity in years
#' @param d Dividend yield, annualized, continuously-compounded
#' @param fn Pricing function name, not in quotes
#' @param f Fully-specified option pricing function, including inputs
#'     which need not be named. For example, you can enter
#'     \code{greeks(bscall(40, 40, .3, .08, .25, 0))}
#' @param ... Pricing function inputs, must be named, may either be a
#'     list or not
#' 
#' @examples
#' s=40; k=40; v=0.30; r=0.08; tt=0.25; d=0;
#' greeks(bscall(s, k, v, r, tt, d))
#' greeks2(bscall, list(s=s, k=k, v=v, r=r, tt=tt, d=d))
#' greeks2(bscall, list(s=s, k=k, v=v, r=r, tt=tt, d=d))[c('Delta', 'Gamma'), ]
#' bsopt(s, k, v, r, tt, d)
#' bsopt(s, c(35, 40, 45), v, r, tt, d)
#' bsopt(s, c(35, 40, 45), v, r, tt, d)[['Call']][c('Delta', 'Gamma'), ]
#'
#' ## plot Greeks for calls and puts for 500 different stock prices
#' ## Unfortunately, this plot will most likely not display in RStudio;
#' ## it will generate a "figure margins too large" error
#' k <- 100; v <- 0.30; r <- 0.08; tt <- 2; d <- 0
#' S <- seq(.5, 250, by=.5)
#' x <- bsopt(S, k, v, r, tt, d)
#' par(mfrow=c(4, 4))  ## create a 4x4 plot
#' for (i in c('Call', 'Put')) {
#'     for (j in rownames(x[[i]])) {  ## loop over greeks
#'         plot(S, x[[i]][j, ], main=paste(i, j), ylab=j, type='l')
#'     }
#' }

#' @export
bsopt <- function(s, k, v, r, tt, d) {
    ## Black-Scholes put and call values.
    xc <- greeks2(bscall, list(s=s, k=k, v=v, r=r, tt=tt, d=d))
    xp <- greeks2(bsput, list(s=s, k=k, v=v, r=r, tt=tt, d=d))
    return(list(Call=xc, Put=xp))
}

#' @export
greeks <- function(f) {
    ## This version uses a standard function call
    args <- match.call()[[2]] ## get f and arguments
    funcname <- as.character(args[[1]])
    args[[1]] <- NULL  ## eliminate function name, leaving only the
                       ## function arguments
    fnames <- names(formals(funcname)) ## arguments to function
    if (sum(names(args)=='') > 0) {
        shared <- intersect(names(args), fnames)
        names(args)[names(args)==''] <- setdiff(fnames, shared)
    } else {
        names(args) <- fnames
    }

    x <<- as.list(args)
    ## Issue: When an argument is a vector, the list representation
    ## stores the values as a language object (an unevaluated
    ## call). In order to extract the values, need to use "eval
    ## x[[i]]".
    for (i in 1:length(x)) x[[i]] <- eval(x[[i]])
    .checkListRecycle(x)
    prem  <-  do.call(funcname, x)
    delta <-  .FirstDer(funcname, 's', x)
    vega  <-  .FirstDer(funcname, 'v', x)/100
    rho   <-  .FirstDer(funcname, 'r', x)/100
    theta <- -.FirstDer(funcname, 'tt', x)/365
    psi   <-  .FirstDer(funcname, 'd', x)/100
    elast <-  x[['s']]*delta/prem
    gamma <-  .SecondDer(funcname, 's', x)
    numcols <- length(prem)
    numrows <- 8
    y <- t(matrix(c(prem,delta,gamma,vega,rho,theta,psi,elast),
                  nrow=numcols,ncol=numrows))
    rownames(y) <- c("Price", "Delta", "Gamma", "Vega", "Rho", "Theta",
                     "Psi", "Elasticity")

    ## In the following, this tests to see if there is variation in
    ## any inputs (is xmaxlength > 1). If so, is there variation in
    ## more than one input (length(maxarg) > 1). The column names are
    ## constructed as appropriate in each case.

    ## are any parameters input as vectors?
    xlength <- lapply(x, length) ## how many of each input?
    xmaxlength <- max(unlist(lapply(x, length))) ## max # of inputs
    arggt1 <- which(xlength > 1)
    if (xmaxlength == 1) {
        colnames(y) <- funcname
    } else {
        ## if we get here, there are multiple inputs with length > 1
        tmp <- NULL
        for (i in arggt1) {
            tmp <- paste(tmp, format(x[[i]], digits=3, trim=TRUE),
                         sep='_')
        }
        colnames(y) <- paste(funcname, tmp, sep='')
    }
    return(y)
}


#' @export
greeks2 <- function(fn, ...) {
    ## Fix handling of inputs with different lengths want to modify
    ## this function so that inputs need not be named
    if (is.list(c(...))) x <- c(...)
    else x <- list(...)
    ## make sure recycling rule will work, stop if not
    .checkListRecycle(x)
    prem  <-  do.call(fn, x)
    delta <-  .FirstDer(fn, 's', x)
    vega  <-  .FirstDer(fn, 'v', x)/100
    rho   <-  .FirstDer(fn, 'r', x)/100
    theta <- -.FirstDer(fn, 'tt', x)/365
    psi   <-  .FirstDer(fn, 'd', x)/100
    elast <-  x[['s']]*delta/prem
    gamma <-  .SecondDer(fn, 1, x)
    numcols <- length(prem)
    numrows <- 8
    y <- t(matrix(c(prem,delta,gamma,vega,rho,theta,psi,elast),
                  nrow=numcols,ncol=numrows))
    rownames(y) <- c("Price", "Delta", "Gamma", "Vega", "Rho", "Theta",
                     "Psi", "Elasticity")
    funcname <- as.character(match.call()[[2]])

    ## The following tests to see if there is variation in any inputs
    ## (is xmaxlength > 1). If so, is there variation in more than one
    ## input (length(maxarg) > 1)? The column names are constructed as
    ## appropriate in each case, showing varying input values by column.

    ## are any parameters input as vectors?
    xlength <- lapply(x, length) ## how many of each input?
    xmaxlength <- max(unlist(lapply(x, length))) ## max # of inputs
    arggt1 <- which(xlength > 1)
    if (xmaxlength == 1) {
        colnames(y) <- funcname
    } else {
        ## if we get here, there are multiple inputs with length > 1
        tmp <- NULL
        for (i in arggt1) {
            tmp <- paste(tmp, format(x[[i]], digits=3, trim=TRUE), sep='_')
        }
        colnames(y) <- paste(funcname, tmp, sep='')
    }
    return(y)
}



.FirstDer <- function(fn, pos, arglist) {
    ## compute first derivative of function fn
    ## arglist must be a list
    epsilon <- 1e-04 
    xup <- xdn <- arglist
    xup[[pos]] <- xup[[pos]] + epsilon
    xdn[[pos]] <- xdn[[pos]] - epsilon
    yup <- do.call(fn, xup)
    ydn <- do.call(fn, xdn)
    return((yup-ydn)/(2*epsilon))
}

.SecondDer <- function(fn, pos, ...) {
    ## this is original
    ##   compute second derivative of function fn
    if (is.list(c(...))) arglist <- c(...)
    else arglist <- list(...)
    epsilon <- 5e-04 
    xup <- xdn <- arglist
    xup[[pos]] <- xup[[pos]] + epsilon
    xdn[[pos]] <- xdn[[pos]] - epsilon
    yup <- .FirstDer(fn, pos, xup)
    ydn <- .FirstDer(fn, pos, xdn)
    return((yup-ydn)/(2*epsilon))
}

.checkListRecycle <- function(x) {
    ## function tests whether list of vectors can work with recylcing
    ## without throwing a warning. We can do this by unlisting the
    ## elements, summing them, and checking for an error. (The summing
    ## will require recycling to work; if it doesn't, there is a
    ## mismatch in the number of entries.)
    tryCatch(
        {tmp <- 0; for (i in seq_along(x)) tmp <- tmp+unlist(x[[i]])},
        warning = function(c) {
            c$message <- paste("Input vector lengths are not",
                               "integer multiples of one another")
            stop(c)
        }
    )
}
