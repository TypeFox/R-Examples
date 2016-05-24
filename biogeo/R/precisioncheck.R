precisioncheck <-
function (dat, x = "x", y = "y", s, e) 
{
    cn <- names(dat)
    f1 <- match(x, cn)
    f2 <- match(y, cn)
    dx <- abs(coord2numeric(dat[, f1]))
    dy <- abs(coord2numeric(dat[, f2]))
    xm <- dx - floor(dx)
    ym <- dy - floor(dy)
    xm <- round(xm, digits = 7)
    ym <- round(ym, digits = 7)
    rs <- seq(s, e, 5)
    cn <- paste("p", rs, "m", sep = "")
    ers <- {
    }
    ersc <- {
    }
    for (j in 1:length(rs)) {
        res <- rs[j]
        mp <- (seq(0, (60 - res), res) + (res/2))/60
        mc <- (seq(0, (60 - res), res))/60
        mp <- round(mp, digits = 7)
        mc <- round(mc, digits = 7)
        er <- rep(2, length(xm))
        erc <- er
        for (i in 1:length(xm)) {
            xp <- (is.element(xm[i], mp)) * 1
            yp <- (is.element(ym[i], mp)) * 1
            xc <- (is.element(xm[i], mc)) * 1
            yc <- (is.element(ym[i], mc)) * 1
            er[i] <- ifelse(yp == 1 & xp == 1, 1, 0)
            erc[i] <- ifelse(yc == 1 & xc == 1, 1, 0)
        }
        ers <- cbind(ers, er)
        ersc <- cbind(ersc, erc)
    }
    ers1 <- as.data.frame(ers)
    ers2 <- as.data.frame(ersc)
    ers3 <- as.data.frame(ers + ersc)
    names(ers1) <- cn
    names(ers2) <- cn
    names(ers3) <- cn
    err <- cbind(dx, dy, ers1)
    err2 <- cbind(dx, dy, ers2)
    preci <- (rowSums(ers3) > 0) * 1
    z <- cbind(dat, ers3, preci)
    return(z)
}
