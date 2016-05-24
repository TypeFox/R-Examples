hc.rma.uni <-
function (object, digits, transf, targs, control, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    x <- object
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    yi <- x$yi
    vi <- x$vi
    k <- length(yi)
    if (k == 1) 
        stop("Stopped because k = 1.")
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
    if (missing(control)) 
        control <- list()
    con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
        verbose = FALSE)
    con[pmatch(names(control), names(con))] <- control
    wi <- 1/vi
    W1 <- sum(wi)
    W2 <- sum(wi^2)/W1
    W3 <- sum(wi^3)/W1
    W4 <- sum(wi^4)/W1
    b <- sum(wi * yi)/W1
    Q <- sum(wi * ((yi - b)^2))
    tau2 <- max(0, (Q - (k - 1))/(W1 - W2))
    vb <- (tau2 * W2 + 1)/W1
    se <- sqrt(vb)
    VR <- 1 + tau2 * W2
    SDR <- sqrt(VR)
    EQ <- function(r) (k - 1) + tau2 * (W1 - W2) + (tau2^2) * 
        ((1/VR^2) * (r^2) - 1/VR) * (W3 - W2^2)
    VQ <- function(r) {
        rsq <- r^2
        recipvr2 <- 1/VR^2
        2 * (k - 1) + 4 * tau2 * (W1 - W2) + 2 * tau2^2 * (W1 * 
            W2 - 2 * W3 + W2^2) + 4 * tau2^2 * (recipvr2 * rsq - 
            1/VR) * (W3 - W2^2) + 4 * tau2^3 * (recipvr2 * rsq - 
            1/VR) * (W4 - 2 * W2 * W3 + W2^3) + 2 * tau2^4 * 
            (recipvr2 - 2 * (1/VR^3) * rsq) * (W3 - W2^2)^2
    }
    scale <- function(r) {
        VQ(r)/EQ(r)
    }
    shape <- function(r) {
        EQ(r)^2/VQ(r)
    }
    finv <- function(f) (W1/W2 - 1) * ((f^2) - 1) + (k - 1)
    Eqn <- function(t) {
        integrand <- function(r) {
            pgamma(finv(r/t), scale = scale(SDR * r), shape = shape(SDR * 
                r)) * dnorm(r)
        }
        integral <- integrate(integrand, lower = t, upper = Inf)$value
        val <- integral - alpha/2
        val
    }
    t0 <- try(uniroot(Eqn, lower = 0, upper = 2, tol = con$tol, 
        maxiter = con$maxiter))
    if (inherits(t0, "try-error")) 
        stop("Error in uniroot().")
    t0 <- t0$root
    u0 <- SDR * t0
    ci.lb <- b - u0 * se
    ci.ub <- b + u0 * se
    b.rma <- x$b
    se.rma <- x$se
    ci.lb.rma <- x$ci.lb
    ci.ub.rma <- x$ci.ub
    if (is.function(transf)) {
        if (is.null(targs)) {
            b <- sapply(b, transf)
            b.rma <- sapply(b.rma, transf)
            se <- NA
            se.rma <- NA
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            ci.lb.rma <- sapply(ci.lb.rma, transf)
            ci.ub.rma <- sapply(ci.ub.rma, transf)
        }
        else {
            b <- sapply(b, transf, targs)
            br.rma <- sapply(b.rma, transf, targs)
            se <- NA
            se.rma <- NA
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            ci.lb.rma <- sapply(ci.lb.rma, transf, targs)
            ci.ub.rma <- sapply(ci.ub.rma, transf, targs)
        }
    }
    res <- list(b = b, se = se, ci.lb = ci.lb, ci.ub = ci.ub, 
        b.rma = b.rma, se.rma = se.rma, ci.lb.rma = ci.lb.rma, 
        ci.ub.rma = ci.ub.rma, method = "DL", method.rma = x$method, 
        tau2 = tau2, tau2.rma = x$tau2, digits = digits)
    class(res) <- c("hc.rma.uni")
    return(res)
}
