

#' @include b_data.R



#' @title Kashp Function
#'
#' @description
#' \code{kashp}, which stands for kappa times arc-sine-hyperbola-power  
#' is the nonlinear transformation of x at the heart 
#' of power hyperbolas, power hyperbolic functions and symmetric Kiener 
#' distributions.
#' \code{dkashp_dx} is its derivative with respect to \code{x}. 
#' \code{ashp} is provided for convenience.
#'
#' @param x a numeric value, vector or matrix.
#' @param k a numeric value or vector, preferably strictly positive.
#'
#' @details
#' \code{ashp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ ashp(x, k) = asinh(x / 2 / k)  }
#' \code{kashp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ kashp(x, k) = k * asinh(x / 2 / k)  }
#' \code{dkashp_dx} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ dkashp_dx(x, k) = 1 / sqrt( x*x/k/k + 4 ) 
#'              = 1 / 2 / cosh( ashp(x, k) ) }
#'
#' If k is a vector, then the use of the function \code{\link[base]{outer}} 
#' is recommanded.
#'
#' The undesired case k=0 returns 0 for kashp and dkashp_dx, 1 for exphp, 
#' -Inf, NaN, +Inf for ashp.
#'
#' @seealso The power hyperbolas and the power hyperbolic functions 
#' \code{\link{exphp}}.
#' 
#' @examples
#' 
#' require(graphics)
#' 
#' ### Example 1
#' x    <- (-3:3)*3 ; names(x) <- x
#' kashp(x, k=2)
#' k    <- c(-2, 0, 1, 2, 3, 5, 10) ; names(k) <- k
#' outer(x, k, kashp)
#' outer(x, k, exphp)
#'
#' ### Example 2
#' xmin   <- -10
#' xd     <- 0.5
#' x      <- seq(xmin, -xmin, xd) ; names(x) <- x
#' k      <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(k) <- k
#' olty   <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd   <- c(1, 1, 2, 2, 3, 4, 2)
#' ocol   <- c(2, 2, 4, 4, 3, 3, 1)
#' op     <- par(mfrow = c(2,2), mgp = c(1.5,0.8,0), mar = c(3,3,2,1))
#' 
#' Tkashp <- ts(cbind(outer(x, k, kashp), "x/2" = x/2), start = xmin, deltat = xd)
#' plot(Tkashp, plot.type = "single", ylim = c(-5, +5), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", xlab = "", 
#'        ylab = "", main = "kashp(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Tkashp), 
#'        inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' Tdkashp <- ts(cbind(outer(x, k, dkashp_dx)), start = xmin, deltat = xd)
#' plot(Tdkashp, plot.type = "single", ylim = c(0, 0.8), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", xlab = "", 
#'        ylab = "", main="dkashp_dx(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Tdkashp), 
#'        inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' Tashp <- ts(cbind(outer(x, k, ashp), "x/2" = x/2), start = xmin, deltat = xd)
#' plot(Tashp, plot.type = "single", ylim = c(-5, +5), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", xlab = "", 
#'        ylab = "", main = "ashp(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Tashp), 
#'        inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End example 2
#' 
#' @export
#' @name kashp
   kashp <- function(x, k = 1) { 
                z <- k * asinh(x / 2 / k) 
                z[which(z == "NaN")] <- 0
                return(z)
               } 
#' @export
#' @rdname kashp
 dkashp_dx <- function(x, k = 1) {
                z <- 1 / sqrt( x*x/k/k + 4 ) 
                z[which(z == "NaN")] <- 0
                return(z)
               } 
#' @export
#' @rdname kashp
    ashp <- function(x, k = 1) { 
                z <- asinh(x / 2 / k) 
                z[which(z == "NaN")] <- 0
                return(z)
               } 







#' @title Power Hyperbolas and Power Hyperbolic Functions
#'
#' @description
#' These functions define the power hyperbola \code{exphp} and the associated 
#' power hyperbolic cosine, sine, tangent, secant, cosecant, cotangent. 
#' They are similar to the traditional hyperbolic functions with term
#' \code{x} receiving a nonlinear transformation via the function 
#' \code{\link{kashp}}.
#'
#' @param x a numeric value, vector or matrix.
#' @param k a numeric value, preferably strictly positive.
#'
#' @details
#' \code{exphp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ exphp(x, k) =  exp( kashp(x, k) )  
#'                          =  exp( k * asinh(x / 2 / k) ) }
#' \code{coshp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ coshp(x, k) = cosh( kashp(x, k) ) }
#' \code{sinhp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ sinhp(x, k) = sinh( kashp(x, k) ) }
#' \code{tanhp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ tanhp(x, k) = tanh( kashp(x, k) ) }
#' \code{sechp} function is defined for x in (-Inf, +Inf) by: 
#'       \deqn{ sechp(x, k) = 1 / coshp(x, k) }
#' \code{cosechp} function is defined for x in (-Inf, 0) U (0, +Inf) by: 
#'       \deqn{ cosechp(x, k) = 1 / sinhp(x, k) }
#' \code{cotanhp} function is defined for x in (-Inf, 0) U (0, +Inf) by: 
#'       \deqn{ cotanhp(x, k) = 1 / tanhp(x, k) }
#'
#' The undesired case k = 0 returns 0 for sinhp and tanhp, 
#' 1 for exphp, coshp and sechp, Inf for cosechp and cotanhp.
#'
#' If k is a vector of length > 1, then the use of the function 
#' \code{\link[base]{outer}} is recommanded.
#'
#' @seealso 
#' The nonlinear transformation \code{\link{kashp}}, the inverse power 
#' hyperbolas and the inverse power hyperbolic functions \code{\link{loghp}}.
#'
#' @examples
#'
#' ### Example 1
#' x  <- (-3:3)*3 
#' exphp(x, k = 4)
#' coshp(x, k = 4)
#' sinhp(x, k = 4) 
#' tanhp(x, k = 4)
#'
#' ### Example 2 outer + plot(exphp, coshp, sinhp, tanhp)
#' xmin  <- -10
#' xd    <- 0.5
#' x     <- seq(xmin, -xmin, xd) ; names(x) <- x
#' k     <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(k) <- k
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 4, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' op    <- par(mfrow = c(2,2), mgp = c(1.5,0.8,0), mar = c(3,3,2,1))
#'
#' ## exphp(x, k)
#' Texphp <- ts(cbind(outer(-x, k, exphp), "exp(-x/2)" = exp(-x/2)), 
#'              start = xmin, deltat = xd)
#' plot(Texphp, plot.type = "single", ylim = c(0,20), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", xlab = "", 
#'        ylab = "", main = "exphp(-x, k)" )
#' legend("topright", title = expression(kappa), legend = colnames(Texphp), 
#'        inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#'
#' ## coshp(x, k)
#' Tcoshp <- ts(cbind(outer(x, k, coshp), "cosh(x/2)" = cosh(x/2)), 
#'              start = xmin, deltat = xd)
#' plot(Tcoshp, plot.type = "single", ylim = c(0,20), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", 
#'        xlab = "", ylab = "", main = "coshp(x, k)" )
#' legend("top", title = expression(kappa), legend = colnames(Tcoshp), 
#'        inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#'
#' ## sinhp(x, k)
#' Tsinhp <- ts(cbind(outer(x, k, sinhp), "sinh(x/2)" = sinh(x/2)), 
#'              start = xmin, deltat=xd)
#' plot(Tsinhp, plot.type = "single", ylim = c(-10,10), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", 
#'        xlab = "", ylab = "", main = "sinhp(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Tsinhp), 
#'        inset = 0.02, lty = olty, lwd=  olwd, col = ocol, cex = 0.7 )
#'
#' ## tanhp(x, k)
#' Ttanhp <- ts(cbind(outer(x, k, tanhp), "tanh(x/2)" = tanh(x/2)), 
#'              start = xmin, deltat = xd)
#' plot(Ttanhp, plot.type = "single", ylim = c(-1,1), 
#'        lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", xlab = "", 
#'        ylab = "", main = "tanhp(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Ttanhp), 
#'        inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End Example 3
#' 
#' @export
#' @name exphp
           exphp <- function(x, k = 1) {  exp( kashp(x, k) ) }
#' @export
#' @rdname exphp
           coshp <- function(x, k = 1) { cosh( kashp(x, k) ) }
#' @export
#' @rdname exphp
           sinhp <- function(x, k = 1) { sinh( kashp(x, k) ) }
#' @export
#' @rdname exphp
           tanhp <- function(x, k = 1) { tanh( kashp(x, k) ) }
#' @export
#' @rdname exphp
           sechp <- function(x, k = 1) { 1 / coshp(x, k) }
#' @export
#' @rdname exphp
         cosechp <- function(x, k = 1) { 1 / sinhp(x, k) }
#' @export
#' @rdname exphp
         cotanhp <- function(x, k = 1) { 1 / tanhp(x, k) }






#' @title Inverse Power Hyperbolas and Inverse Power Hyperbolic Functions
#'
#' @description
#' The inverse power hyperbolas and the inverse power hyperbolic functions: 
#' arc-cosine-hp, arc-sine-hp, arc-tangent-hp, 
#' arc-secant-hp, arc-cosecant-hp and arc-cotangent-hp.
#'
#' @param x a numeric value, vector or matrix.
#' @param k a numeric value, preferably strictly positive.
#'
#' @details
#' \code{loghp} function is defined on (0, +Inf) by: 
#'       \deqn{  logshp(x, k) = 2 * k * sinh( log(x) / k) }
#' \code{acoshp} function is defined on [1, +Inf) by: 
#'       \deqn{ acoshp(x, k) = 2 * k * sinh( acosh(x) / k) }
#' \code{asinhp} function is defined on (-Inf, +Inf) by: 
#'       \deqn{ asinhp(x, k) = 2 * k * sinh( asinh(x) / k) }
#' \code{atanhp} function is defined on (-1, +1) by: 
#'       \deqn{ atanhp(x, k) = 2 * k * sinh( atanh(x) / k) }
#' \code{asechp} function is defined on (0, +1] by: 
#'       \deqn{ asechp(x, k) = 2 * k * sinh( acosh(1/x) / k) }
#' \code{acosechp} function is defined on (-Inf, 0) U (0, +Inf) by: 
#'     \deqn{ acosechp(x, k) = 2 * k * sinh( asinh(1/x) / k) }
#' \code{acotanhp} function is defined on (-Inf, -1) U (1, +Inf) by: 
#'     \deqn{ acotanhp(x, k) = 2 * k * sinh( atanh(1/x) / k) }
#'
#' If k is a vector of length > 1, then the use of the function 
#' \code{\link[base]{outer}} is recommanded.
#'
#'
#' @seealso The power hyperbolic functions \code{\link{exphp}}.
#' @examples
#' 
#' ### Example 1 (acoshp, asinhp, atanhp)
#'  loghp( c(ppoints(10), 1, 1/rev(ppoints(10))), k = 2)
#' acoshp( 1:10, k = 2)
#' asinhp( -5:5, k = 2)
#' atanhp( seq(-1, 1, by = 0.1), k = 2)
#' asechp( ppoints(20), k = 2)
#' acosechp( -5:5, k = 2)
#' acotanhp( c( -1/ppoints(10), 1/rev(ppoints(10))), k = 2)
#'
#' x  <- (-3:3)*3 
#'  loghp(exphp(x, k = 4), k = 4)
#' acoshp(coshp(x, k = 4), k = 4)
#' asinhp(sinhp(x, k = 4), k = 4)
#' atanhp(tanhp(x, k = 4), k = 4)
#'
#' 
#' ### Example 2 (loghp, acoshp, asinhp, atanhp)
#' k     <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(k) <- k
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 4, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' op    <- par(mfrow = c(2, 2), mgp = c(1.5, 0.8, 0), mar = c(3, 3, 2, 1))
#' 
## loghp(x, k)
#' xld     <- 0.05
#' xl      <- seq(0.05, 20, xld) ; names(xl) <- xl
#' Tlcoshp <- ts(cbind(outer(xl, k, loghp), "2*log(x)" = 2*log(xl)), 
#'               start = xl[1], deltat = xld)
#' plot(Tlcoshp, plot.type = "single", xlim = c(0,20), ylim = c(-5,15), 
#'      lty = olty, lwd = olwd, col = ocol, xaxs = "i", yaxs = "i", 
#'      xlab="", ylab = "", main = "loghp(x, k)" )
#' legend("bottomright", title = expression(kappa), legend = colnames(Tlcoshp), 
#'      inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#'
#' ## acoshp(x, k)
#' xcd     <- 0.5
#' xc      <- seq(1, 20, xcd) ; names(xc) <- xc
#' Tacoshp <- ts(cbind(outer(xc, k, acoshp), "2*acosh(x)" = 2*acosh(xc)), 
#'               start = xc[1], deltat = xcd)
#' plot(Tacoshp, plot.type = "single", ylim = c(0,15), lty = olty, lwd = olwd, col = ocol,
#'         xaxs = "i", yaxs = "i", xlab = "", ylab = "", main = "acoshp(x, k)" )
#' legend("bottomright", title = expression(kappa), legend = colnames(Tacoshp), 
#'         inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' ## asinhp(x, k)
#' xsd     <- 0.5
#' xs      <- seq(-10, 10, xsd) ; names(xs) <- xs
#' Tasinhp <- ts(cbind(outer(xs, k, asinhp), "2*asinh(x)" = 2*asinh(xs)), 
#'               start = xs[1], deltat = xsd)
#' plot(Tasinhp, plot.type = "single", ylim = c(-10,10), lty = olty, lwd = olwd, col = ocol,
#'         xaxs = "i", yaxs = "i", xlab = "", ylab = "", main = "asinhp(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Tasinhp), 
#'         inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' ## atanhp(x, k)
#' xtd     <- 0.01
#' xt      <- seq(-1, 1, xtd) ; names(xt) <- xt
#' Tatanhp <- ts(cbind(outer(xt, k, atanhp), "2*atanh(x)" = 2*atanh(xt)), 
#'               start = xt[1], deltat = xtd)
#' plot(Tatanhp, plot.type = "single", ylim = c(-10,10), lty = olty, lwd = olwd, col = ocol,
#'         xaxs = "i", yaxs = "i", xlab = "", ylab = "", main = "atanhp(x, k)" )
#' legend("topleft", title = expression(kappa), legend = colnames(Tatanhp), 
#'         inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End Example 2
#' 
#'
#' @export
#' @name loghp
   loghp <- function(x, k = 1) { 2 * k * sinh( log(x) / k) }
#' @export
#' @rdname loghp
  acoshp <- function(x, k = 1) { 2 * k * sinh( acosh(x) / k) }
#' @export
#' @rdname loghp
  asinhp <- function(x, k = 1) { 2 * k * sinh( asinh(x) / k) }
#' @export
#' @rdname loghp
  atanhp <- function(x, k = 1) { 2 * k * sinh( atanh(x) / k) }
#' @export
#' @rdname loghp
  asechp <- function(x, k = 1) { 2 * k * sinh( acosh(1/x) / k) }
#' @export
#' @rdname loghp
acosechp <- function(x, k = 1) { 2 * k * sinh( asinh(1/x) / k) }
#' @export
#' @rdname loghp
acotanhp <- function(x, k = 1) { 2 * k * sinh( atanh(1/x) / k) }



