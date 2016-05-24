SSgampt <-
structure(function (input, ks, A) 
{
    .value <- ks + (A/input)
    .actualArgs <- as.list(match.call()[c("ks", "A")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 2), list(NULL, c("ks", 
            "A")))
        .grad[, "ks"] <- 1
        .grad[, "A"] <- 1/input
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .value
}, initial = function (mCall, data, LHS) 
{
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    if (nrow(xy) < 3) {
        stop("Too few distinct input values to fit a Philip model")
    }
    ndistinct <- nrow(xy)
    fc <- mean(xy[(ndistinct - 2):ndistinct, ][["y"]])
    lfirst <- xy[1:(ndistinct/4), ]
    pars2 <- coef(lm(exp(y) ~ sqrt(x), data = lfirst))
    S <- (-1/pars2[2])
    value <- c(ks = fc, A = S)
    names(value) <- mCall[c("ks", "A")]
    value
}, pnames = c("ks", "A"), class = "selfStart")
