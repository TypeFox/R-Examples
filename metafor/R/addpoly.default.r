addpoly.default <-
function (x, vi, sei, ci.lb, ci.ub, rows = -1, level = 95, digits = 2, 
    annotate = TRUE, mlab, transf, atransf, targs, efac = 1, 
    col, border, cex, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(transf)) 
        transf <- FALSE
    if (missing(atransf)) 
        atransf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (missing(mlab)) 
        mlab <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(col)) 
        col <- "black"
    if (missing(border)) 
        border <- "black"
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    yi <- x
    if (hasArg(ci.lb) && hasArg(ci.ub)) {
        if (length(ci.lb) != length(ci.ub)) 
            stop("Length of ci.lb and ci.ub do not match.")
        if (missing(vi) && missing(sei)) {
            vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
        }
        else {
            if (missing(vi)) 
                vi <- sei^2
        }
        if (length(ci.lb) != length(vi)) 
            stop("Length of vi (or sei) does not match length of (ci.lb, ci.ub) pairs.")
    }
    else {
        if (missing(vi)) {
            if (missing(sei)) {
                stop("Must specify either vi, sei, or (ci.lb, ci.ub) pairs.")
            }
            else {
                vi <- sei^2
                ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
                ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
            }
        }
        else {
            ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
            ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
        }
    }
    k <- length(yi)
    if (is.null(rows)) {
        rows <- -1:(-k)
    }
    else {
        if (length(rows) == 1L) {
            rows <- rows:(rows - k + 1)
        }
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'rows' argument.")
    yivi.na <- is.na(yi) | is.na(vi)
    if (any(yivi.na)) {
        not.na <- !yivi.na
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ci.lb <- ci.lb[not.na]
            ci.ub <- ci.ub[not.na]
            mlab <- mlab[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq_len(length(rows.na))) {
                rows.new[rows <= rows.na[j]] <- rows.new[rows <= 
                  rows.na[j]] + 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    yi.ut <- yi
    ci.lb.ut <- ci.lb
    ci.ub.ut <- ci.ub
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    cex.adj <- min(1, 20/height)
    xlim <- par.usr[1:2]
    if (is.null(cex)) 
        cex <- par("cex") * cex.adj
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- round(cbind(sapply(yi.ut, atransf), 
                  sapply(ci.lb.ut, atransf), sapply(ci.ub.ut, 
                    atransf)), digits)
            }
            else {
                annotext <- round(cbind(sapply(yi.ut, atransf, 
                  targs), sapply(ci.lb.ut, atransf, targs), sapply(ci.ub.ut, 
                  atransf, targs)), digits)
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            annotext <- round(cbind(yi, ci.lb, ci.ub), digits)
        }
        annotext <- matrix(apply(annotext, 2, format, nsmall = digits), 
            ncol = 3)
        annotext <- cbind(annotext[, 1], " [ ", annotext[, 2], 
            " , ", annotext[, 3], " ]")
        annotext <- apply(annotext, 1, paste, collapse = "")
        text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, 
            ...)
    }
    for (i in seq_len(k)) {
        polygon(x = c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y = c(rows[i], 
            rows[i] + (height/100) * cex * efac, rows[i], rows[i] - 
                (height/100) * cex * efac), col = col, border = border, 
            ...)
        if (!is.null(mlab)) 
            text(xlim[1], rows[i], mlab[i], pos = 4, cex = cex, 
                ...)
    }
}
