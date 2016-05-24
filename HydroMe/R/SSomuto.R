SSomuto <-
structure(function (input, Ths1, alp1, Ths2, alp2) 
{
    .expr1 <- exp(alp1)
    .expr4 <- exp(((-.expr1) * input))
    .expr7 <- exp(alp2)
    .expr10 <- exp(((-.expr7) * input))
    .value <- (Ths1 * .expr4) + (Ths2 * .expr10)
    .actualArgs <- as.list(match.call()[c("Ths1", "alp1", "Ths2", 
        "alp2")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 4), list(NULL, c("Ths1", 
            "alp1", "Ths2", "alp2")))
        .grad[, "Ths1"] <- .expr4
        .grad[, "alp1"] <- -(Ths1 * (.expr4 * (.expr1 * input)))
        .grad[, "Ths2"] <- .expr10
        .grad[, "alp2"] <- -(Ths2 * (.expr10 * (.expr7 * input)))
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .value
}, initial = function (mCall, data, LHS) 
{
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    if (nrow(xy) < 5) {
        stop("Too few distinct input values to fit a biexponential")
    }
    ndistinct <- nrow(xy)
    nlast <- max(3, round(ndistinct/2))
    dfirst1 <- xy[1, ][["y"]][1]
    dlast1 <- xy[nrow(xy), ][["y"]][1]
    B <- ifelse(((xy[1, ][["y"]] - xy[ndistinct, ][["y"]]) < 
        (xy[1, ][["y"]])/2), 0, xy[ndistinct, ][["y"]])
    dlast <- xy[(ndistinct + 1 - nlast):ndistinct, ]
    pars2 <- coef(lm(logb(y) ~ x, data = dlast))
    lrc2 <- logb(abs(pars2[2]))
    xy[["res"]] <- xy[["y"]] - exp(pars2[1]) * exp(-exp(lrc2) * 
        xy[["x"]])
    dfirst <- xy[1:(ndistinct - nlast), ]
    pars1 <- coef(lm(logb(abs(res)) ~ x, data = dfirst))
    lrc1 <- logb(abs(pars1[2]))
    pars <- coef(nls(y ~ cbind(exp(-exp(lrc1) * x), exp(-exp(lrc2) * 
        x)), data = xy, start = list(lrc1 = lrc1, lrc2 = lrc2), 
        algorithm = "plinear"))
    value <- c(Ths1 = pars[3], alp1 = pars[1], Ths2 = pars[4], 
        alp2 = pars[2])
    names(value) <- mCall[c("Ths1", "alp1", "Ths2", "alp2")]
    value
}, pnames = c("Ths1", "alp1", "Ths2", "alp2"), class = "selfStart")
