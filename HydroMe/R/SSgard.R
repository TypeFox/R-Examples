SSgard <-
structure(function (input, Thr, Ths, alp, scal) 
{
    .expr1 <- Ths - Thr
    .expr2 <- input^scal
    .expr4 <- 1 + (alp * .expr2)
    .expr5 <- .expr4^-1
    .expr9 <- .expr4^-2
    .value <- Thr + (.expr1 * .expr5)
    .actualArgs <- as.list(match.call()[c("Thr", "Ths", "alp", 
        "scal")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 4), list(NULL, c("Thr", 
            "Ths", "alp", "scal")))
        .grad[, "Thr"] <- 1 - .expr5
        .grad[, "Ths"] <- .expr5
        .grad[, "alp"] <- .expr1 * (-1 * (.expr2 * .expr9))
        .grad[, "scal"] <- .expr1 * (-1 * ((alp * (.expr2 * (log(input)))) * 
            .expr9))
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .value
}, initial = function (mCall, data, LHS) 
{
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    if (nrow(xy) < 5) {
        stop("Too few distinct input values to fit a van Genuchten model")
    }
    ndistinct <- nrow(xy)
    nlast <- max(3, round(ndistinct/2))
    Thr <- ifelse(((xy[1, ][["y"]] - xy[ndistinct, ][["y"]]) < 
        (xy[1, ][["y"]])/2), 0, xy[ndistinct, ][["y"]])
    Ths <- xy[1:(ndistinct - nlast), ][["y"]][1]
    dmid <- xy[(ndistinct/2 - 1):(ndistinct/2 + 2), ]
    pars2 <- coef(lm(y ~ log(x), data = dmid))
    scal <- exp(pars2[2])
    alp <- 1/xy[(ndistinct/2 - 1):(ndistinct/2 + 2), ][["x"]][1]
    value <- c(Thr = Thr, Ths = Ths, alp = alp, scal = scal)
    names(value) <- mCall[c("Thr", "Ths", "alp", "scal")]
    value
}, pnames = c("Thr", "Ths", "alp", "scal"), class = "selfStart")
