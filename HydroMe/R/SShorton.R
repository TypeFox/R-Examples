SShorton <-
structure(function (input, fc, f0, lrk) 
{
    .expr1 <- f0 - fc
    .expr2 <- exp(lrk)
    .expr5 <- exp(((-.expr2) * input))
    .value <- fc + (.expr1 * .expr5)
    .actualArgs <- as.list(match.call()[c("fc", "f0", "lrk")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 3), list(NULL, c("fc", 
            "f0", "lrk")))
        .grad[, "fc"] <- 1 - .expr5
        .grad[, "f0"] <- .expr5
        .grad[, "lrk"] <- -(.expr1 * (.expr5 * (.expr2 * input)))
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .value
}, initial = function (mCall, data, LHS) 
{
    xy <- sortedXyData(mCall[["input"]], LHS, data)
    if (nrow(xy) < 3) {
        stop("Too few distinct input values to fit a asymptotic regression model")
    }
    if (nrow(xy) > 3) {
        xy$ydiff <- abs(xy$y - NLSstRtAsymptote(xy))
        xy <- data.frame(xy)
        lrc <- logb(-coef(lm(logb(ydiff) ~ x, data = xy))[2])
        names(lrc) <- NULL
        pars <- coef(nls(y ~ cbind(1 - exp(-exp(lrc) * x), exp(-exp(lrc) * 
            x)), data = xy, start = list(lrc = lrc), algorithm = "plinear"))
    }
    else {
        ydiff <- diff(xy$y)
        if (prod(ydiff) <= 0) {
            stop("Can't fit an asymptotic regression model to these data")
        }
        avg.resp <- xy$y
        frac <- (avg.resp[3] - avg.resp[1])/(avg.resp[2] - avg.resp[1])
        xunique <- unique(xy$x)
        xdiff <- diff(xunique)
        if (xdiff[1] == xdiff[2]) {
            expmRd <- frac - 1
            rc <- -logb(expmRd)/xdiff[1]
            lrc <- logb(rc)
            expmRx1 <- exp(-rc * xunique[1])
            bma <- ydiff[1]/(expmRx1 * (expmRd - 1))
            Asym <- avg.resp[1] - bma * expmRx1
            pars <- c(lrk = lrc, fc = Asym, f0 = bma + Asym)
        }
        else {
            stop("Too few observations to fit an asymptotic regression model")
        }
    }
    names(pars) <- NULL
    val <- list(pars[2], pars[3], pars[1])
    names(val) <- mCall[c("fc", "f0", "lrk")]
    val
}, pnames = c("fc", "f0", "lrk"), class = "selfStart")
