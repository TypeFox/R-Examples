#' @title The Elliptic Distribution
#'
#' @description
#' Density, distribution function, quantile function, 
#' and random generation for the univariate elliptic distribution.
#'
#' @param p numeric vector of probabilities.
#' @param x numeric vector of quantiles.
#' @param q numeric vector of quantiles.
#' @param n number of observations.
#' @param object an object of ecd class. 
#'               To achieve high performance for \code{qec} and \code{rec}, 
#'               it should be created with \code{with.quantile=TRUE}.
#' @param debug logical, whether to print debug message, default is \code{FALSE}.
#'
#' @return \code{dec} gives the density,
#' \code{pec} gives the distribution function, 
#' \code{qec} gives the quantile function, 
#' \code{rec} generates random deviates.
#'
#' @keywords distribution
#'
#' @include ecd-class.R
#'
#' @export dec 
#' @export pec 
#' @export qec 
#' @export rec
#'
#' @importFrom stats predict runif
#' 
#' @author Stephen H. Lihn
#'
#' @examples
#' d <- ecd(with.quantile=TRUE)
#' x <- seq(-20, 20, by=5)
#' dec(x,d)
#' pec(x,d)
#' p <- c(0.0001, 0.001, 0.01, 0.99, 0.999, 0.9999)
#' qec(p,d)
#' rec(100,d)
### <======================================================================>
"dec" <- function(x, object=ecd()) {
    ecd.pdf(object, x)
}

#' @rdname dec
"pec" <- function(q, object = ecd()) {

    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    N <- 4
    mean <- object@stats$mean
    dv <- object@stats$stdev
    xt1 <- mean - N*dv # left tail
    xt2 <- mean + N*dv # right tail
    
    # (a) switch piece.wise off for tail region 
    # (b) switch to ccdf for q > mean
    f1 <- function(q) ecd.cdf(object, q, piece.wise=TRUE)
    f2 <- function(q) { 1-ecd.ccdf(object, q, piece.wise=TRUE) }
    f3 <- function(q) ecd.cdf(object, q, piece.wise=FALSE)
    f4 <- function(q) { 1-ecd.ccdf(object, q, piece.wise=FALSE) }
    
    Nq <- length(q)
    p <- ecd.mpnum(object, rep(NaN, Nq))
    insert <- function(i,v) { p[i] <<- v }
    
    i1 <- which(q <= mean & q > xt1)
    i2 <- which(q > mean & q < xt2)
    i3 <- which(q <= xt1)
    i4 <- which(q >= xt2)
    
    if (length(i1)>0) mapply(insert, i1, f1(q[i1]))
    if (length(i2)>0) mapply(insert, i2, f2(q[i2]))
    if (length(i3)>0) mapply(insert, i3, f3(q[i3]))
    if (length(i4)>0) mapply(insert, i4, f4(q[i4]))
    ecd.mpnum(object, p)
}

#' @rdname dec
"qec" <- function(p, object = ecd(with.quantile=TRUE), debug=FALSE) {

    if (! ecd.has_quantile(object)) { 
        object <- quantilize(object, show.warning=TRUE) 
    }
    dq <- object@quantile

    p <- ecd.mp2f(p) # lm can't handle mpfr, so revert to numeric
    
    in_range <- (p>=0 & p<=1)
    if (! all(in_range)) {
        stop("Probability input out of range (0,1)!")
    }
    
    between_seg <- function(idx, c.from, c.to, q) {
        if (q >= c.from & q <= c.to) idx else 0
    }

    # translate q to segment index
    get_idx <- function(q) {
        max(mapply(between_seg, 
                   seq(1, dq@N_seg), 
                   dq@cseg.from, dq@cseg.to, 
                   MoreArgs = list(q=q), 
                   SIMPLIFY = TRUE))
    }
    
    predict_quantile <- function(q) {
        idx <- get_idx(q)
        if (debug) print(paste("idx", idx, "of", dq@N_seg, "for q=", q))
        
        if (idx > 0) {
            if (debug) {
                print(paste("xseg from", dq@xseg.from[idx], " to ", dq@xseg.to[idx]))
                print(paste("cseg from", dq@cseg.from[idx], " to ", dq@cseg.to[idx]))
            }
            fit <- dq@cdf.fit[[idx]]
            if (debug) print(fit)           
            xp <- predict(fit, data.frame(cdf=c(q)))
        } else {
            # tail region
            if (q <= dq@cseg.min & q > 0) {
                # left tail
                if (debug) {
                    print(paste("left tail, log cdf for q=", q))
                    print(dq@fit.left)
                }
                xp <- predict(dq@fit.left, data.frame(lcdf=c(log(q))))
            } else if (q >= dq@cseg.max & q < 1) {
                # right tail
                if (debug) {
                    print(paste("right tail, log ccdf for 1-q, q=", q))
                    print(dq@fit.right)
                }
                xp <- predict(dq@fit.right, data.frame(lccdf=c(log(1-q))))
            } 
            else if (q==0) { # edge case
                if (debug) print(paste("edge case for q=", q))
                xp <- -Inf 
            } else if (q==1) { 
                if (debug) print(paste("edge case for q=", q))
                xp <- Inf
            } else {
                stop(paste("Failed to locate quantile region. q=",q))
            } 
        }
        xp
    }
    
    # quantile specified, give me x
    if (debug & length(p)==1) {
        return(ecd.mpnum(object, predict_quantile(p)))
    }
    
    rs <- simplify2array(parallel::mclapply(p, predict_quantile))
    ecd.mpnum(object, rs)
}

#' @rdname dec
"rec" <- function(n, object = ecd(with.quantile=TRUE)) {
    qec(runif(n), object)
}
### <---------------------------------------------------------------------->
