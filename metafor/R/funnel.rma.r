funnel.rma <-
function (x, yaxis = "sei", xlim, ylim, xlab, ylab, steps = 5, 
    at, atransf, targs, digits, level = x$level, addtau2 = FALSE, 
    type = "rstandard", back = "lightgray", shade = "white", 
    hlines = "white", refline, pch = 19, pch.fill = 21, ci.res = 1000, 
    ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (is.element("robust.rma", class(x))) 
        stop("Function not applicable to objects of class \"robust.rma\".")
    na.act <- getOption("na.action")
    yaxis <- match.arg(yaxis, c("sei", "vi", "seinv", "vinv", 
        "ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi"))
    type <- match.arg(type, c("rstandard", "rstudent"))
    if (missing(atransf)) 
        atransf <- FALSE
    atransf.char <- deparse(substitute(atransf))
    if (missing(ylab)) {
        if (yaxis == "sei") 
            ylab <- "Standard Error"
        if (yaxis == "vi") 
            ylab <- "Variance"
        if (yaxis == "seinv") 
            ylab <- "Inverse Standard Error"
        if (yaxis == "vinv") 
            ylab <- "Inverse Variance"
        if (yaxis == "ni") 
            ylab <- "Sample Size"
        if (yaxis == "ninv") 
            ylab <- "Inverse Sample Size"
        if (yaxis == "sqrtni") 
            ylab <- "Square-Root Sample Size"
        if (yaxis == "sqrtninv") 
            ylab <- "Inverse Square-Root Sample Size"
        if (yaxis == "lni") 
            ylab <- "Log Sample Size"
        if (yaxis == "wi") 
            ylab <- "Weight (in %)"
    }
    if (missing(at)) 
        at <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", 
        "lni"))) {
        if (is.null(x$ni)) 
            stop("No sample size information stored in model object.")
        if (anyNA(x$ni)) 
            warning("Sample size information stored in model object \n  contains NAs. Not all studies will be plotted.")
    }
    if (missing(digits)) {
        if (yaxis == "sei") 
            digits <- c(2, 3)
        if (yaxis == "vi") 
            digits <- c(2, 3)
        if (yaxis == "seinv") 
            digits <- c(2, 3)
        if (yaxis == "vinv") 
            digits <- c(2, 3)
        if (yaxis == "ni") 
            digits <- c(2, 0)
        if (yaxis == "ninv") 
            digits <- c(2, 3)
        if (yaxis == "sqrtni") 
            digits <- c(2, 3)
        if (yaxis == "sqrtninv") 
            digits <- c(2, 3)
        if (yaxis == "lni") 
            digits <- c(2, 3)
        if (yaxis == "wi") 
            digits <- c(2, 2)
    }
    else {
        if (length(digits) == 1L) 
            digits <- c(digits, digits)
    }
    if (x$int.only) {
        if (missing(refline)) 
            refline <- x$b
        if (is.element("rma.mv", class(x))) {
            addtau2 <- FALSE
            tau2 <- 0
        }
        tau2 <- ifelse(addtau2, x$tau2, 0)
        yi <- x$yi
        vi <- x$vi
        ni <- x$ni
        sei <- sqrt(vi)
        slab <- x$slab[x$not.na]
        if (missing(xlab)) 
            xlab <- .setlab(x$measure, transf.char = "FALSE", 
                atransf.char, gentype = 1)
    }
    else {
        if (missing(refline)) 
            refline <- 0
        tau2 <- 0
        options(na.action = "na.pass")
        if (type == "rstandard") {
            res <- rstandard(x)
        }
        else {
            res <- rstudent(x)
        }
        options(na.action = na.act)
        not.na <- !is.na(res$resid)
        yi <- res$resid[not.na]
        sei <- res$se[not.na]
        ni <- x$ni.f[not.na]
        vi <- sei^2
        slab <- x$slab[not.na]
        if (missing(xlab)) 
            xlab <- "Residual Value"
    }
    if (yaxis == "wi") {
        options(na.action = "na.omit")
        weights <- weights(x)
        options(na.action = na.act)
    }
    if (missing(ylim)) {
        if (yaxis == "sei") 
            ylim <- c(max(sei), 0)
        if (yaxis == "vi") 
            ylim <- c(max(vi), 0)
        if (yaxis == "seinv") 
            ylim <- c(min(1/sei), max(1/sei))
        if (yaxis == "vinv") 
            ylim <- c(min(1/vi), max(1/vi))
        if (yaxis == "ni") 
            ylim <- c(min(ni, na.rm = TRUE), max(ni, na.rm = TRUE))
        if (yaxis == "ninv") 
            ylim <- c(max(1/ni, na.rm = TRUE), min(1/ni, na.rm = TRUE))
        if (yaxis == "sqrtni") 
            ylim <- c(min(sqrt(ni), na.rm = TRUE), max(sqrt(ni), 
                na.rm = TRUE))
        if (yaxis == "sqrtninv") 
            ylim <- c(max(1/sqrt(ni), na.rm = TRUE), min(1/sqrt(ni), 
                na.rm = TRUE))
        if (yaxis == "lni") 
            ylim <- c(min(log(ni), na.rm = TRUE), max(log(ni), 
                na.rm = TRUE))
        if (yaxis == "wi") 
            ylim <- c(min(weights), max(weights))
    }
    else {
        if (is.element(yaxis, c("sei", "vi", "ninv", "sqrtninv"))) 
            ylim <- c(max(ylim), min(ylim))
        if (is.element(yaxis, c("seinv", "vinv", "ni", "sqrtni", 
            "lni", "wi"))) 
            ylim <- c(min(ylim), max(ylim))
        if (is.element(yaxis, c("sei", "vi", "ni", "ninv", "sqrtni", 
            "sqrtninv", "lni"))) {
            if (ylim[1] < 0 || ylim[2] < 0) 
                stop("Both limits for the y axis must be >= 0.")
        }
        if (is.element(yaxis, c("seinv", "vinv"))) {
            if (ylim[1] <= 0 || ylim[2] <= 0) 
                stop("Both limits for the y axis must be > 0.")
        }
        if (is.element(yaxis, c("wi"))) {
            if (ylim[1] < 0 || ylim[2] < 0) 
                stop("Both limits for the y axis must be >= 0.")
        }
    }
    if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {
        alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
        alpha.min <- min(alpha)
        avals <- length(alpha)
        if (yaxis == "sei") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1]^2 + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1]^2 + tau2)
        }
        if (yaxis == "vi") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1] + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(ylim[1] + tau2)
        }
        if (yaxis == "seinv") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1]^2 + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1]^2 + tau2)
        }
        if (yaxis == "vinv") {
            x.lb.bot <- refline - qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1] + tau2)
            x.ub.bot <- refline + qnorm(alpha.min/2, lower.tail = FALSE) * 
                sqrt(1/ylim[1] + tau2)
        }
        if (missing(xlim)) {
            xlim <- c(min(x.lb.bot, min(yi)), max(x.ub.bot, max(yi)))
            rxlim <- xlim[2] - xlim[1]
            xlim[1] <- xlim[1] - (rxlim * 0.1)
            xlim[2] <- xlim[2] + (rxlim * 0.1)
        }
        else {
            xlim <- sort(xlim)
        }
    }
    if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", 
        "lni", "wi"))) {
        if (missing(xlim)) {
            xlim <- c(min(yi), max(yi))
            rxlim <- xlim[2] - xlim[1]
            xlim[1] <- xlim[1] - (rxlim * 0.1)
            xlim[2] <- xlim[2] + (rxlim * 0.1)
        }
        else {
            xlim <- sort(xlim)
        }
    }
    if (!is.null(at)) {
        xlim[1] <- min(c(xlim[1], at), na.rm = TRUE)
        xlim[2] <- max(c(xlim[2], at), na.rm = TRUE)
    }
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
        xaxt = "n", yaxt = "n", bty = "n", ...)
    par.usr <- par("usr")
    rect(par.usr[1], par.usr[3], par.usr[2], par.usr[4], col = back, 
        border = NA, ...)
    axis(side = 2, at = seq(from = ylim[1], to = ylim[2], length.out = steps), 
        labels = formatC(seq(from = ylim[1], to = ylim[2], length.out = steps), 
            digits = digits[2], format = "f"), ...)
    abline(h = seq(from = ylim[1], to = ylim[2], length.out = steps), 
        col = hlines, ...)
    if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {
        if (yaxis == "sei") {
            rylim <- ylim[1] - ylim[2]
            ylim[1] <- ylim[1] + (rylim * 0.1)
            ylim[2] <- max(0, ylim[2] - (rylim * 0.1))
        }
        if (yaxis == "vi") {
            rylim <- ylim[1] - ylim[2]
            ylim[1] <- ylim[1] + (rylim * 0.1)
            ylim[2] <- max(0, ylim[2] - (rylim * 0.1))
        }
        if (yaxis == "seinv") {
            rylim <- ylim[2] - ylim[1]
            ylim[2] <- ylim[2] + (rylim * 0.1)
        }
        if (yaxis == "vinv") {
            rylim <- ylim[2] - ylim[1]
            ylim[2] <- ylim[2] + (rylim * 0.1)
        }
        yi.vals <- seq(from = ylim[1], to = ylim[2], length.out = ci.res)
        if (yaxis == "sei") 
            vi.vals <- yi.vals^2
        if (yaxis == "vi") 
            vi.vals <- yi.vals
        if (yaxis == "seinv") 
            vi.vals <- 1/yi.vals^2
        if (yaxis == "vinv") 
            vi.vals <- 1/yi.vals
        for (m in avals:1) {
            ci.left <- refline - qnorm(alpha[m]/2, lower.tail = FALSE) * 
                sqrt(vi.vals + tau2)
            ci.right <- refline + qnorm(alpha[m]/2, lower.tail = FALSE) * 
                sqrt(vi.vals + tau2)
            polygon(c(ci.left, ci.right[ci.res:1]), c(yi.vals, 
                yi.vals[ci.res:1]), border = NA, col = shade[m], 
                ...)
            lines(ci.left, yi.vals, lty = "dotted", ...)
            lines(ci.right, yi.vals, lty = "dotted", ...)
        }
    }
    if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) 
        segments(refline, ylim[1], refline, ylim[2], ...)
    if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", 
        "lni", "wi"))) 
        abline(v = refline, ...)
    xaxis.vals <- yi
    if (yaxis == "sei") 
        yaxis.vals <- sei
    if (yaxis == "vi") 
        yaxis.vals <- vi
    if (yaxis == "seinv") 
        yaxis.vals <- 1/sei
    if (yaxis == "vinv") 
        yaxis.vals <- 1/vi
    if (yaxis == "ni") 
        yaxis.vals <- ni
    if (yaxis == "ninv") 
        yaxis.vals <- 1/ni
    if (yaxis == "sqrtni") 
        yaxis.vals <- sqrt(ni)
    if (yaxis == "sqrtninv") 
        yaxis.vals <- 1/sqrt(ni)
    if (yaxis == "lni") 
        yaxis.vals <- log(ni)
    if (yaxis == "wi") 
        yaxis.vals <- weights
    points(xaxis.vals, yaxis.vals, pch = pch, ...)
    if (is.element("rma.uni.trimfill", class(x))) 
        points(xaxis.vals[x$fill], yaxis.vals[x$fill], pch = pch.fill, 
            col = "black", bg = "white", ...)
    box(bty = "l")
    if (is.null(at)) {
        at <- axTicks(side = 1)
    }
    else {
        at <- at[at > par("usr")[1]]
        at <- at[at < par("usr")[2]]
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[1], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[1], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[1], format = "f")
    }
    axis(side = 1, at = at, labels = at.lab, ...)
    sav <- data.frame(x = xaxis.vals, y = yaxis.vals, slab = slab)
    sav$fill <- x$fill
    invisible(sav)
}
