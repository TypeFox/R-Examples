SSfredlund <-
structure(function (x, thr, ths, alp, nscal, mscal) 
{
    .expr1 <- ths - thr
    .expr2 <- x/alp
    .expr3 <- .expr2^nscal
    .expr4 <- 2.7183 + .expr3
    .expr5 <- log(.expr4)
    .expr6 <- .expr5^mscal
    .expr9 <- 1/.expr6
    .expr12 <- .expr5^(mscal - 1)
    .expr23 <- .expr6^2
    .value <- thr + .expr1/.expr6
    .grad <- array(0, c(length(.value), 5L), list(NULL, c("thr", 
        "ths", "alp", "nscal", "mscal")))
    .grad[, "thr"] <- 1 - .expr9
    .grad[, "ths"] <- .expr9
    .grad[, "alp"] <- .expr1 * (.expr12 * (mscal * (.expr2^(nscal - 
        1) * (nscal * (x/alp^2))/.expr4)))/.expr23
    .grad[, "nscal"] <- -(.expr1 * (.expr12 * (mscal * (.expr3 * 
        log(.expr2)/.expr4)))/.expr23)
    .grad[, "mscal"] <- -(.expr1 * (.expr6 * log(.expr5))/.expr23)
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 4) {
      stop("Too few distinct x values to fit a Fredlund model")
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
    pars <- as.numeric(c(thr=Thr1,ths = Ths1, alp = 1/alp1, nscal=scal1, mscal=m1))
    setNames(c(pars[1], pars[2], pars[3], pars[4], pars[5]), mCall[c("thr","ths","alp","nscal","mscal")])

}, pnames = c("thr", "ths", "alp", "nscal", "mscal"), class = "selfStart")
