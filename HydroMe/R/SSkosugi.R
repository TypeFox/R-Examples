SSkosugi <-
structure(function (x, thr, ths, alp, nscal) 
{
    .expr2 <- (ths - thr) * 2
    .expr3 <- x * alp
    .expr4 <- log(.expr3)
    .expr6 <- sqrt(2)
    .expr7 <- .expr4/nscal * .expr6
    .expr8 <- pnorm(.expr7)
    .expr11 <- 2 * .expr8
    .expr13 <- dnorm(.expr7)
    .value <- thr + .expr2 * .expr8
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("thr", 
        "ths", "alp", "nscal")))
    .grad[, "thr"] <- 1 - .expr11
    .grad[, "ths"] <- .expr11
    .grad[, "alp"] <- .expr2 * (.expr13 * (x/.expr3/nscal * .expr6))
    .grad[, "nscal"] <- -(.expr2 * (.expr13 * (.expr4/nscal^2 * 
        .expr6)))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 4) {
      stop("Too few distinct x values to fit a Kosugi model")
    }
    ndistinct <- nrow(xy)
    nlast <- max(3, round(ndistinct/2))
    dfirst <- xy[1, ][["y"]][1]
    dlast <- xy[nrow(xy), ][["y"]][1]
    Thr1 <- ifelse(((xy[1, ][["y"]] - xy[ndistinct, ][["y"]]) <(xy[1, ][["y"]])/2), 0.01, xy[ndistinct, ][["y"]])
    Ths1 <- dfirst
    dmid <- xy[(ndistinct/2 - 2):(ndistinct/2 + 1), ]
    pars2 <- coef(lm(y ~ log(x), data = dmid))
    ymid <- xy[1:max(3, round(nrow(xy)/2)), ][["y"]][max(3, round(nrow(xy)/2))]
    ax <- (ymid - pars2[1])/pars2[2]
    slopep <- pars2[2]/(dfirst - dlast)
    m1 <- ifelse(abs(slopep) < 1, (1 - exp(-0.8 * (abs(slopep)))),(1 - (0.5755/(abs(slopep))) + (0.1/(abs(slopep))^2) +(0.025/(abs(slopep))^3)))
    scal1 <- 1/(1 - m1)
    alp1 <- (((2^(1/m1)) - 1)^(1 - m1))/exp(ax)
    pars <- as.numeric(c(thr=Thr1,ths = Ths1, alp = 1/alp1, nscal=scal1))
    setNames(c(pars[1], pars[2], pars[3], pars[4]), mCall[c("thr","ths","alp","nscal")])

}, pnames = c("thr", "ths", "alp", "nscal"), class = "selfStart")
