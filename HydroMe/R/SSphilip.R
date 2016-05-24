SSphilip <-
structure(function (input, fc, S) 
{
    .expr2 <- input^-0.5
    .value <- fc + ((0.5 * S) * .expr2)
    .actualArgs <- as.list(match.call()[c("fc", "S")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 2), list(NULL, c("fc", 
            "S")))
        .grad[, "fc"] <- 1
        .grad[, "S"] <- 0.5 * .expr2
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
    value <- c(fc = fc, S = S)
    names(value) <- mCall[c("fc", "S")]
    value
}, pnames = c("fc", "S"), class = "selfStart")
