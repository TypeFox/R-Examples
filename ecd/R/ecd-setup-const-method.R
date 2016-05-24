#' Integration preprocessor for an ecd object
#'
#' This is an internal helper function for ecd constructor
#' Its main function is to determine \code{const}, \code{const_left_x}, and \code{const_right_x} during object construction.
#'
#' @param object An object of ecd class
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return \code{list(const, const_left_x, const_right_x)}
#'
#' @keywords constructor
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @examples
#' ecd.toString(ecd(-1,1, sigma=0.1))
#'
### <======================================================================>
ecd.setup_const <- function(object, verbose=FALSE) {
    # This is a sub-system under ecd constructor.
    
    # It figures out the constant for the PDF
    
    pdf2 <- function(x) {
	    y <- solve(object, x)
	    exp(y)
	}
    
    R <- object@R
    degree <- floor(ecd.mp2f(object@theta/pi*180))
    sigma <- ecd.mp2f(object@sigma)
    step <- 2*sigma
    
    # special performance tweak for small angles
    if (object@use.mpfr & degree < 10 & R < 8) {
        step <- sigma
    }

    tol <- 10* .Machine$double.eps^0.25 * sigma
    if (R > 4) {
        tol <- tol/sigma * ecd.mp2f(ecd.estimate_const(object))
    }
    
    c1 <- ecd.mp2f(object@mu - step)
    c2 <- ecd.mp2f(object@mu + step)
    max_retry <- 10
    
    total_iter <- 0
    
    # ---------------------------------------------------------------------------------
    # left wing
    if (verbose) print(paste(Sys.time(), "ecd.setup_const: left wing, deg=", degree, "step=", step, "tol=", tol))

    sd <- object@R^(1/3)*object@sigma
    d1 <- object@mu - 4*sd
    d2 <- object@mu + 4*sd
    sd4_special <- FALSE
    if (object@use.mpfr) {
        if (object@R > 1000 & tol < .Machine$double.eps) { # 1e-16
            sd4_special <- TRUE
        }
        if (object@gamma==0 & object@alpha <=-100) {
            sd4_special <- TRUE
        }
    }    
    retry <- max_retry
    p1 <- list(message="Empty")
	repeat {
        if (sd4_special) { 
            if (verbose) print(paste(Sys.time(), "ecd.setup_const: split 4*sd"))
            p1  <- ecd.integrate(object, pdf2, -Inf, d1, abs.tol=tol, show.warning=FALSE)
            p1a <- ecd.integrate(object, pdf2,   d1, c1, abs.tol=tol, show.warning=FALSE)
            p1 <- .ecd.intg_merge(p1, p1a)
            
        } else {
            p1 <- ecd.integrate(object, pdf2, -Inf, c1, abs.tol=tol, show.warning=FALSE)
        }
        if (p1$message == "OK") break
        c1 <- c1 - step
        retry <- retry - 1
        if ( retry < 0 ) {
            stop(paste("Failed to integrate left wing of pdf:",
                       "c1=",c1, "msg=", p1$message))
        }
        if (verbose) print(paste(Sys.time(), "ecd.setup_const: retry", max_retry-retry, "c1", c1))
    }
    total_iter <- total_iter + (max_retry-retry)
    if (p1$message != "OK") .ecd.setup_const_warn(object, "left wing const is not OK.", p1)
    if (!(p1$value>=0)) {
        stop(paste("Failed to integrate left wing of pdf:",
                   "c1=",c1, "p1=", p1$value))
    }
    if (verbose) print(paste(Sys.time(), "ecd.setup_const: left wing, abs.error=", 
                             ecd.mp2f(p1$abs.error)))
    
    # ---------------------------------------------------------------------------------
    # right wing
    if (verbose) print(paste(Sys.time(), "ecd.setup_const: right wing"))
    retry <- max_retry
    p2 <- list(message="Empty")
	repeat {
        if (sd4_special) { 
            if (verbose) print(paste(Sys.time(), "ecd.setup_const: split 4*sd"))
            p2  <- ecd.integrate(object, pdf2, c2,  d2, abs.tol=tol, show.warning=FALSE)
            p2a <- ecd.integrate(object, pdf2, d2, Inf, abs.tol=tol, show.warning=FALSE)
            p2 <- .ecd.intg_merge(p2, p2a)
            
        } else {
            p2 <- ecd.integrate(object, pdf2, c2, Inf, abs.tol=tol, show.warning=FALSE)
        }
        if (p2$message == "OK") break
        c2 <- c2 + step
        retry <- retry - 1
        if ( retry < 0 ) {
            stop(paste("Failed to integrate right wing of pdf:",
                       "c2=",c2, "msg=", p2$message))
        }
        if (verbose) print(paste(Sys.time(), "ecd.setup_const: retry", max_retry-retry, "c2", c2))
    }
    total_iter <- total_iter + (max_retry-retry)
    if (p2$message != "OK") .ecd.setup_const_warn(object, "right wing const is not OK.", p2)
    if (!(p2$value>=0)) {
        stop(paste("Failed to integrate right wing of pdf:",
                   "c2=",c2, "p2=", p2$value))
    }    
    if (verbose) print(paste(Sys.time(), "ecd.setup_const: right wing, abs.error=", 
                             ecd.mp2f(p2$abs.error)))
    
    # ---------------------------------------------------------------------------------
    # center
    if (verbose) print(paste(Sys.time(), "ecd.setup_const: center"))
    
	p3 <- ecd.integrate(object, pdf2, c1, object@mu, abs.tol=tol, show.warning=FALSE)
    if (p3$message != "OK") .ecd.setup_const_warn(object, "center left const is not OK.", p3)
    if (!(p3$value>=0)) {
        stop(paste("Failed to integrate center left of pdf:",
                   "c1=", c1, "c2=", object@mu, "p3=", p3$value))
    }

    p4 <- ecd.integrate(object, pdf2, object@mu, c2, abs.tol=tol, show.warning=FALSE)
    if (p4$message != "OK") .ecd.setup_const_warn(object, "center right const is not OK.", p4)
    if (!(p4$value>=0)) {
        stop(paste("Failed to integrate center right of pdf:",
                   "c1=", object@mu, "c2=", c2, "p3=", p4$value))
    }
    if (verbose) print(paste(Sys.time(), "ecd.setup_const: center, abs.error=", 
                             ecd.mp2f(p3$abs.error+p4$abs.error)))
    
    # ---------------------------------------------------------------------------------
    # sum it up
    p_total <- p1$value + p2$value + p3$value + p4$value
    if (verbose) {
        print(paste(Sys.time(), "p1=", ecd.mp2f(p1$value), "p2=", ecd.mp2f(p2$value),
                    "p3=", ecd.mp2f(p3$value), "p4=", ecd.mp2f(p4$value)
        ))
        print(paste(Sys.time(), "ecd.setup_const: done",
                    "c1=", ecd.mp2f(c1), "c2=", ecd.mp2f(c2), "const=", ecd.mp2f(p_total)
        ))
    }
    list(const = unname(p_total),
         const_left_x = unname(c1),
         const_right_x = unname(c2),
         total_iter = total_iter
        )
}
### <---------------------------------------------------------------------->
.ecd.setup_const_warn <- function(object, message, p) {
    warning(paste("WARN: setup_const-", message, "(", p$message, ")", "from ecd:", ecd.toString(object)))
}
