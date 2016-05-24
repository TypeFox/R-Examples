coef.glmpath.cr <-
function (object, s, weight = NULL,  
    eps = .Machine$double.eps, ...) 
{
    mode <- "step"
    type <- "coefficients"
    b <- object$b.corrector
    std.b <- scale(b[, -1], FALSE, 1/object$sdx)
    if (!is.null(object$nopenalty.subset)) 
        std.b <- std.b[, -object$nopenalty.subset, drop = FALSE]
    k <- nrow(b)
    steps <- seq(k)
    if (missing(s)) {
        s <- steps[object$new.A]
        if (mode != "step") {
            warning("no s argument; mode switched to step")
            mode <- "step"
        }
    }
    sb <- switch(mode, step = {
        if (any(s < 1) | any(s > k)) 
            stop("Argument s out of range")
        steps
    })
    sfrac <- (s - sb[1])/(sb[k] - sb[1])
    sb <- (sb - sb[1])/(sb[k] - sb[1])
    usb <- unique(sb)
    useq <- match(usb, sb)
    sb <- sb[useq]
    b <- b[useq, ]
    coord <- approx(sb, seq(sb), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newb <- ((sb[right] - sfrac) * b[left, , drop = FALSE] + 
        (sfrac - sb[left]) * b[right, , drop = FALSE])/(sb[right] - 
        sb[left])
    newb[left == right, ] <- b[left[left == right], ]
    fit <- newb
    dimnames(fit) <- list(s, object$xnames)
    attr(fit, "s") <- s
    attr(fit, "fraction") <- sfrac
    attr(fit, "mode") <- mode
    fit
}

