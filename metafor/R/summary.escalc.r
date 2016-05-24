summary.escalc <-
function (object, out.names = c("sei", "zi", "ci.lb", "ci.ub"), 
    var.names, H0 = 0, append = TRUE, replace = TRUE, level = 95, 
    digits, transf, ...) 
{
    if (!is.element("escalc", class(object))) 
        stop("Argument 'object' must be an object of class \"escalc\".")
    x <- object
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    if (length(out.names) != 4) 
        stop("Argument out.names must be of length 4.")
    if (any(out.names != make.names(out.names, unique = TRUE))) {
        out.names <- make.names(out.names, unique = TRUE)
        warning(paste0("Argument 'out.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: out.names = c('", 
            out.names[1], "', '", out.names[2], "', '", out.names[3], 
            "', '", out.names[4], "')."))
    }
    if (missing(transf)) 
        transf <- FALSE
    if (missing(var.names)) {
        if (!is.null(attr(x, "yi.names"))) {
            yi.name <- attr(x, "yi.names")[1]
        }
        else {
            if (!is.element("yi", names(x))) 
                stop("Cannot determine name of the 'yi' variable.")
            yi.name <- "yi"
        }
        if (!is.null(attr(x, "vi.names"))) {
            vi.name <- attr(x, "vi.names")[1]
        }
        else {
            if (!is.element("vi", names(x))) 
                stop("Cannot determine name of the 'vi' variable.")
            vi.name <- "vi"
        }
    }
    else {
        if (length(var.names) != 2) 
            stop("Argument 'var.names' must be of length 2.")
        if (any(var.names != make.names(var.names, unique = TRUE))) {
            var.names <- make.names(var.names, unique = TRUE)
            warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", 
                var.names[1], "', '", var.names[2], "')."))
        }
        yi.name <- var.names[1]
        vi.name <- var.names[2]
    }
    yi <- x[[yi.name]]
    vi <- x[[vi.name]]
    if (is.null(yi) || is.null(vi)) 
        stop(paste0("Cannot find variables '", yi.name, "' and/or '", 
            vi.name, "' in the data frame."))
    k <- length(yi)
    if (length(H0) == 1) 
        H0 <- rep(H0, k)
    sei <- sqrt(vi)
    zi <- (yi - H0)/sei
    if (is.function(transf)) {
        ci.lb <- mapply(transf, yi - crit * sei, ...)
        ci.ub <- mapply(transf, yi + crit * sei, ...)
        yi <- mapply(transf, yi, ...)
        attr(x, "transf") <- TRUE
    }
    else {
        ci.lb <- yi - crit * sei
        ci.ub <- yi + crit * sei
        attr(x, "transf") <- FALSE
    }
    x[[yi.name]] <- yi
    x[[vi.name]] <- vi
    if (append) {
        dat <- data.frame(x)
        if (replace) {
            dat[[out.names[1]]] <- sei
            dat[[out.names[2]]] <- zi
            dat[[out.names[3]]] <- ci.lb
            dat[[out.names[4]]] <- ci.ub
        }
        else {
            if (is.element(out.names[1], names(dat))) {
                is.na.sei <- is.na(dat[[out.names[1]]])
                dat[[out.names[1]]][is.na.sei] <- sei[is.na.sei]
            }
            else {
                dat[[out.names[1]]] <- sei
            }
            if (is.element(out.names[2], names(dat))) {
                is.na.zi <- is.na(dat[[out.names[2]]])
                dat[[out.names[2]]][is.na.zi] <- zi[is.na.zi]
            }
            else {
                dat[[out.names[2]]] <- zi
            }
            if (is.element(out.names[3], names(dat))) {
                is.na.ci.lb <- is.na(dat[[out.names[3]]])
                dat[[out.names[3]]][is.na.ci.lb] <- ci.lb[is.na.ci.lb]
            }
            else {
                dat[[out.names[3]]] <- ci.lb
            }
            if (is.element(out.names[4], names(dat))) {
                is.na.ci.ub <- is.na(dat[[out.names[4]]])
                dat[[out.names[4]]][is.na.ci.ub] <- ci.ub[is.na.ci.ub]
            }
            else {
                dat[[out.names[4]]] <- ci.ub
            }
        }
    }
    else {
        dat <- data.frame(yi, vi, sei, zi, ci.lb, ci.ub)
        names(dat) <- c(yi.name, vi.name, out.names)
    }
    if (!missing(digits)) {
        attr(dat, "digits") <- digits
    }
    else {
        attr(dat, "digits") <- attr(x, "digits")
    }
    if (is.null(attr(dat, "digits"))) 
        attr(dat, "digits") <- 4
    if (!missing(var.names)) {
        attr(dat, "yi.names") <- unique(c(var.names[1], attr(object, 
            "yi.names")))
    }
    else {
        attr(dat, "yi.names") <- unique(c(yi.name, attr(object, 
            "yi.names")))
    }
    if (!missing(var.names)) {
        attr(dat, "vi.names") <- unique(c(var.names[2], attr(object, 
            "vi.names")))
    }
    else {
        attr(dat, "vi.names") <- unique(c(vi.name, attr(object, 
            "vi.names")))
    }
    attr(dat, "sei.names") <- unique(c(out.names[1], attr(object, 
        "sei.names")))
    attr(dat, "zi.names") <- unique(c(out.names[2], attr(object, 
        "zi.names")))
    attr(dat, "ci.lb.names") <- unique(c(out.names[3], attr(object, 
        "ci.lb.names")))
    attr(dat, "ci.ub.names") <- unique(c(out.names[4], attr(object, 
        "ci.ub.names")))
    class(dat) <- c("escalc", "data.frame")
    return(dat)
}
