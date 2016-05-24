`ICEalice` <-
function (ICEw) 
{
    if (missing(ICEw) || !inherits(ICEw, "ICEwedge")) 
        stop("The first argument to ICEalice must be an ICEwedge object.")
    lambda <- ICEw$lambda
    unit <- ICEw$unit
    ia <- as.vector(ICEw$axys[, 1])
    n <- length(ia)
    a11 <- c(45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 
        120, 127.5, 135)
    wtp <- round(lambda * tan((a11 - 45)/57.29578), digits = 3)
    wtp[13] <- Inf
    wta <- round(lambda/tan((a11 - 45)/57.29578), digits = 3)
    wta[1] <- Inf
    acc <- as.matrix(cbind(a11, wtp, rep(0, 13), wta, rep(0, 
        13)))
    dimnames(acc) <- list(rep(1:13), c("ICEangle", "WTP", "VAGR", 
        "WTA", "ALICE"))
    neqcl <- swqcl <- 0    
    for (i in 1:n) {
        for (j in 1:13) {
            if (ia[i] <= acc[j, 1] && ia[i] >= acc[j, 1] - 180) 
                acc[j, 3] <- acc[j, 3] + 1
            if (ia[i] <= acc[j, 1] && ia[i] >= -acc[j, 1]) 
                acc[j, 5] <- acc[j, 5] + 1
        }
        if (ia[i] > 45 && ia[i] <= 135) 
            neqcl <- neqcl + 1
        if (ia[i] < -45 && ia[i] >= -135) 
            swqcl <- swqcl + 1     
    }
    acc[, 3] <- acc[, 3]/n
    acc[, 5] <- acc[, 5]/n
    neqcl <- neqcl/n
    swqcl <- swqcl/n
    qcl <- round(100 * as.vector(cbind(acc[1,5], neqcl, swqcl, 1-acc[13,5])), digits = 1)
    ICEaccol <- list(lambda = lambda, unit = unit, ia = ia, acc = acc, qcl = qcl)
    class(ICEaccol) <- "ICEalice"
    ICEaccol
}
`ICEangle` <-
function (t1, t2) 
{
    if (length(t1) == 0) 
        t1 <- 0
    if (length(t2) == 0) 
        t2 <- 0
    u <- t1 - t2
    v <- t1 + t2
    ia <- 90
    if (abs(u) > 1e-05 * abs(v)) 
        ia <- atan(abs(v)/u) * 57.29578
    if (ia < 0) 
        ia <- ia + 180
    if (v < 0) 
        ia <- -ia
    ia
}
`ICEcolor` <-
function (ICEw, lfact = 1, beta = 1, gamma = 3 + 2 * sqrt(2)) 
{
    if (missing(ICEw) || !inherits(ICEw, "ICEwedge")) 
        stop("The first argument to ICEcolor must be an existing ICEwedge object.")
    if (lfact <= 0) 
        stop("The lfact argument to ICEcolor must be strictly positive.")
    if (beta <= 0) 
        stop("The beta argument to ICEcolor must be strictly positive.")
    if (gamma <= 0) 
        stop("The gamma = eta*beta argument to ICEcolor must be strictly positive.")
    lambda <- lfact * ICEw$lambda
    t1 <- ICEw$t1
    unit <- ICEw$unit
    conf <- ICEw$conf
    axys <- ICEw$axys
    xmax <- ICEw$xmax
    ymax <- ICEw$ymax
    if (lfact != 1) {
        if (unit == "cost") {
            t1[1] <- t1[1] * lfact
            axys[, 2] <- axys[, 2] * lfact
            xmax <- xmax * lfact
        }
        else {
            t1[2] <- t1[2]/lfact
            axys[, 3] <- axys[, 3]/lfact
            ymax <- ymax/lfact
        }
    }
    r <- as.vector(sqrt(axys[, 2]^2 + axys[, 3]^2))
    dif <- (axys[, 2] - axys[, 3])/sqrt(2)
    abcos <- abs(dif)/r
    pref <- sign(dif) * r^beta * abcos^gamma
    ICEclwol <- list(lambda = lambda, beta = beta, gamma = gamma, 
        unit = unit, axys = axys, conf = conf, pref = pref, xmax = xmax, 
        ymax = ymax, jlo = ICEw$jlo, kup = ICEw$kup)
    class(ICEclwol) <- "ICEcolor"
    ICEclwol
}
`ICEd2m` <-
function (d, rownew, rowstd) 
{
    e0 <- mean(d[rowstd, 2])
    e1 <- mean(d[rownew, 2])
    c0 <- mean(d[rowstd, 3])
    c1 <- mean(d[rownew, 3])
    c(e1 - e0, c1 - c0)
}
`ICEepmap` <-
function (lambda = 1, beta = 1, gamma = 3 + 2 * sqrt(2)) 
{
    if (lambda <= 0) 
        stop("The Lambda argument to ICEepmap must be strictly positive.")
    if (beta <= 0) 
        stop("The Beta argument to ICEepmap must be strictly positive.")
    if (gamma <= 0) 
        stop("The Gamma = eta*beta argument to ICEepmap must be strictly positive.")
    ICEepmol <- list(lambda = lambda, beta = beta, gamma = gamma)
    class(ICEepmol) <- "ICEepmap"
    ICEepmol
}
`ICEomega` <-
function (lambda = 1, beta = 1, eta = 3 + 2 * sqrt(2)) 
{
    if (lambda <= 0) 
        stop("The Lambda argument to ICEomega must be strictly positive.")
    if (beta <= 0) 
        stop("The Beta argument to ICEomega must be strictly positive.")
    if (eta <= 0) 
        stop("The Eta = gamma/beta ratio argument to ICEomega must be strictly positive.")
    gamma = eta * beta
    ICEepmol <- list(lambda = lambda, beta = beta, gamma = gamma)
    class(ICEepmol) <- "ICEepmap"
    ICEepmol
}
`ICErunif` <-
function (n1, n2) 
{
    m <- floor(n2 + 1 - n1)
    idx <- as.integer(floor(m * runif(m)) + floor(n1))
    idx
}
`ICEscale` <-
function (df, trtm, xeffe, ycost, lambda = 1, unit = "cost") 
{
    if (missing(df) || !inherits(df, "data.frame")) 
        stop("The first argument to ICEscale must be an existing Data Frame.")
    if (missing(trtm)) 
        stop("The Second argument to ICEscale must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(df)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(df[, trtm])) != 2) 
        stop("Treatment factor must assume exactly two different levels.")
    if (missing(xeffe)) 
        stop("The Third argument to ICEscale must name the Treatment Effectiveness variable.")
    xeffe <- deparse(substitute(xeffe))
    if (!is.element(xeffe, dimnames(df)[[2]])) 
        stop("Effectiveness measure must be an existing Data Frame variable.")
    if (missing(ycost)) 
        stop("The Fourth argument to ICEscale must name the Treatment Cost variable.")
    ycost <- deparse(substitute(ycost))
    if (!is.element(ycost, dimnames(df)[[2]])) 
        stop("Cost measure must be an existing Data Frame variable.")
    if (lambda <= 0) 
        stop("The lambda argument to ICEscale must be strictly positive.")
    if (unit != "effe") 
        unit <- "cost"
    effcst <- na.omit(df[, c(trtm, xeffe, ycost)])
    if (unit != "cost") 
        effcst[, 3] <- effcst[, 3]/lambda
    else effcst[, 2] <- effcst[, 2] * lambda
    names(effcst) <- c("trtm", "effe", "cost")
    effcst <- effcst[do.call(order, effcst), ]
    idx <- table(as.numeric(effcst$trtm))
    rowstd <- 1:idx[1]
    rownew <- 1:idx[2] + idx[1]
    nstd <- length(rowstd)
    nnew <- length(rownew)
    t1 <- c(mean(effcst[rownew, 2]) - mean(effcst[rowstd, 2]), 
        mean(effcst[rownew, 3]) - mean(effcst[rowstd, 3]))
    s1 <- c(sqrt(var(effcst[rownew, 2])/nnew + var(effcst[rowstd, 
        2])/nstd), sqrt(var(effcst[rownew, 3])/nnew + var(effcst[rowstd, 
        3])/nstd))
    ICEsclol <- list(trtm = trtm, xeffe = xeffe, ycost = ycost, 
        effcst = effcst, lambda = lambda, unit = unit, t1 = t1, 
        s1 = s1)
    class(ICEsclol) <- "ICEscale"
    ICEsclol
}
`ICEuncrt` <-
function (df, trtm, xeffe, ycost, lambda = 1, unit = "cost", 
    R = 25000, seed = 0) 
{
    if (missing(df) || !inherits(df, "data.frame")) 
        stop("The first argument to ICEuncrt must be an existing Data Frame.")
    if (lambda <= 0) 
        stop("The lambda argument to ICEuncrt must be strictly positive.")
    if (unit != "effe") 
        unit <- "cost"
    ICEuncol <- list(df = deparse(substitute(df)), lambda = lambda, unit = unit, R = R)
    if (missing(trtm)) 
        stop("The Second argument to ICEuncrt must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(df)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(df[, trtm])) != 2) 
        stop("Treatment factor must assume exactly two different levels.")
    if (missing(xeffe)) 
        stop("The Third argument to ICEuncrt must name the Treatment Effectiveness variable.")
    xeffe <- deparse(substitute(xeffe))
    if (!is.element(xeffe, dimnames(df)[[2]])) 
        stop("Effectiveness measure must be an existing Data Frame variable.")
    if (missing(ycost)) 
        stop("The Fourth argument to ICEuncrt must name the Treatment Cost variable.")
    ycost <- deparse(substitute(ycost))
    if (!is.element(ycost, dimnames(df)[[2]])) 
        stop("Cost measure must be an existing Data Frame variable.")
    effcst <- na.omit(df[, c(trtm, xeffe, ycost)])
    names(effcst) <- c("trtm", "effe", "cost")
    effcst <- effcst[do.call(order, effcst), ]
    if (unit != "cost") 
        effcst[, 3] <- effcst[, 3]/lambda
    else effcst[, 2] <- effcst[, 2] * lambda
    if (R > 25000) 
        R <- 25000
    else if (R < 50) 
        R <- 50
    idx <- table(as.numeric(effcst$trtm))
    rowstd <- 1:idx[1]
    rownew <- 1:idx[2] + idx[1]
    nstd <- length(rowstd)
    nnew <- nstd + length(rownew)
    t1 <- ICEd2m(effcst, rownew, rowstd)
    t <- matrix(rep(0, R * 2), nrow = R, ncol = 2)
    if (seed == 0)
        seed <- 1 + floor(25000 * runif(1))
    set.seed(seed)
    for (i in 1:R) {
        rowstd <- ICErunif(1, nstd)
        rownew <- ICErunif(nstd + 1, nnew)
        t[i, ] <- ICEd2m(effcst, rownew, rowstd)
    }
    ICEuncol <- c(ICEuncol, list(trtm = trtm, xeffe = xeffe, 
        ycost = ycost, effcst = effcst, t1 = t1, t = t, seed = seed))
    class(ICEuncol) <- "ICEuncrt"
    ICEuncol
}
`ICEwedge` <-
function (ICEu, lfact = 1, conf = 0.95) 
{
    if (missing(ICEu) || !inherits(ICEu, "ICEuncrt")) 
        stop("The first argument to ICEwedge must be an ICEuncrt object.")
    if (lfact < 0)
        stop("The lfact argument to ICEwedge must be non-negative.")
    if (conf < 0.5 || conf > 0.99) 
        stop("Wedge Confidence Level must be within [0.50, 0.99].")
    unit <- ICEu$unit
    R <- ICEu$R
    t <- ICEu$t
    t1 <- ICEu$t1
    ab <- "alibi"
    if (lfact == 0) {
        lfact <- max(c(abs(max(t[, 2])), abs(min(t[, 2])))) /
                 max(c(abs(max(t[, 1])), abs(min(t[, 1]))))
        ab <- "alias"                
        }
    lambda <- lfact * ICEu$lambda
    if (lfact != 1) {
        if (unit == "cost") {
            t1[1] <- t1[1] * lfact
            t[, 1] <- t[, 1] * lfact
        }
        else {
            t1[2] <- t1[2]/lfact
            t[, 2] <- t[, 2]/lfact
        }
    }
    ICEwdgol <- list(ICEinp = deparse(substitute(ICEu)), lambda = lambda, 
        lfact = lfact, unit = unit, conf = conf)
    ia1 <- ICEangle(t1[1], t1[2])
    ia <- rep(90, R)
    for (i in 1:R) ia[i] <- ICEangle(t[i, 1], t[i, 2])
    axys <- data.frame(cbind(ia, t, rep(0, R)))
    axys <- axys[do.call(order, axys), ]
    if (ia1 < 0) {
        if (ia1 < axys[1, 1]) 
            center <- 1
        else {
            for (j in 1:(R - 1)) {
                i <- j + 1
                if (axys[j, 1] <= ia1 && ia1 < axys[i, 1]) {
                  center <- j
                  break
                }
            }
        }
    }
    else {
        if (ia1 > axys[R, 1]) 
            center <- R
        else {
            for (j in (R - 1):1) {
                i <- j + 1
                if (axys[j, 1] < ia1 && ia1 <= axys[i, 1]) {
                  center <- j
                  break
                }
            }
        }
    }
    j <- center - floor(R * conf/2)
    if (j < 1) 
        j <- j + R
    k <- j + floor(R * conf)
    if (k > R) 
        k <- k - R
    subangle <- axys[k, 1] - axys[j, 1]
    if (subangle < 0) 
        subangle <- subangle + 360
    if (j < k) 
        axys[j:k, 4] <- 1
    else {
        axys[1:k, 4] <- 1
        axys[j:R, 4] <- 1
    }
    xmax <- max(c(abs(max(axys[, 2])), abs(min(axys[, 2]))))
    ymax <- max(c(abs(max(axys[, 3])), abs(min(axys[, 3]))))
    ICEwdgol <- c(ICEwdgol, list(R = R, axys = axys, t1 = t1, 
        ia1 = ia1, center = center, jlo = j, kup = k, subangle = subangle, 
        xmax = xmax, ymax = ymax, ab = ab))
    class(ICEwdgol) <- "ICEwedge"
    ICEwdgol
}
`plot.ICEalice` <-
function (x, ...) 
{
    if (missing(x) || !inherits(x, "ICEalice")) 
        stop("The first argument to plot.ICEalice must be an ICEalice object.")
    a7 <- c(45, 60, 75, 90, 105, 120, 135)
    x1 <- x$acc[, 1]
    x2 <- x$acc[, 2]
    x2[13] <- 4 * x2[12]
    y3 <- x$acc[, 3]
    y5 <- x$acc[, 5]
    matplot(x2, y3, type = "n", xlim = c(0, 3 * x2[12]), xlab = "ICE Willingness To Pay", 
        ylim = c(0, 1), ylab = "Acceptability", main = "VAGR Acceptability Curve")
    matlines(x2, y3)
    opar <- par(ask = dev.interactive(orNone = TRUE))
    matplot(x1, y5, type = "n", xlim = c(45, 135), xlab = "ICE Angle (theta)", 
        ylim = c(0, 1), ylab = "Acceptability", main = "ALICE Curve", 
        axes = FALSE)
    axis(1, a7)
    axis(2)
    box()
    matlines(x1, y5)
    abline(v = 90)
    par(opar)
}
`plot.ICEcolor` <-
function (x, alibi = FALSE, ...) 
{
    if (missing(x) || !inherits(x, "ICEcolor")) 
        stop("The first argument to plot.ICEcolor must be an ICEcolor object.")
    cv <- rainbow(12, start = 0, end = 0.33)
    pmax <- max(c(abs(max(x$pref, na.rm = TRUE)), abs(min(x$pref, 
        na.rm = TRUE))))
    eta <- x$gamma/x$beta
    if (alibi == FALSE) {
        plot(x$axys[x$axys[, 4] == 1, 2], x$axys[x$axys[, 4] == 1, 
            3], main = "ICE Alias Confidence Wedge with Preference Colors", 
            xlab = "Effectiveness Difference", ylab = "Cost Difference", 
            sub = paste("lambda =", round(x$lambda, digits = 3), ", beta =", 
                round(x$beta, digits = 3), ", gamma =", round(x$gamma, 
                digits = 3), ", eta =", round(eta, digits = 3)), 
            xlim = c(-1 * x$xmax, x$xmax), ylim = c(-1 * x$ymax, x$ymax),
            pch = 20, bg = "white", col = cv[round(5.5 * (x$pref[x$axys[, 
                4] == 1] + pmax)/pmax) + 1])
        points(x$axys[x$axys[, 4] == 0, 2], x$axys[x$axys[, 4] == 0, 3],
            col = "black", pch = 20)
        par(lty = 1)
        abline(h = 0, v = 0)
        mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$kup, 2]^2 + 
            x$axys[x$kup, 3]^2)
        xray <- c(0, x$axys[x$kup, 2]) * mfac
        yray <- c(0, x$axys[x$kup, 3]) * mfac
        lines(xray, yray)
        mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$jlo, 2]^2 + 
            x$axys[x$jlo, 3]^2)
        xray <- c(0, x$axys[x$jlo, 2]) * mfac
        yray <- c(0, x$axys[x$jlo, 3]) * mfac
        lines(xray, yray)
        par(lty = 3)
        abline(c(0, 1))
        }
    else {
        amax = max(x$xmax, x$ymax)
        plot(x$axys[x$axys[, 4] == 1, 2], x$axys[x$axys[, 4] == 1, 
            3], main = "ICE Alibi Confidence Wedge with Preference Colors", 
            xlab = "Effectiveness Difference", ylab = "Cost Difference", 
            sub = paste("lambda =", round(x$lambda, digits = 3), ", beta =", 
                round(x$beta, digits = 3), ", gamma =", round(x$gamma, 
                digits = 3), ", eta =", round(eta, digits = 3)), 
            xlim = c(-amax, amax), asp = 1, pch = 20, bg = "white",
            col = cv[round(5.5 * (x$pref[x$axys[, 4] == 1] + pmax)/pmax) + 1])
        points(x$axys[x$axys[, 4] == 0, 2], x$axys[x$axys[, 4] == 0, 3],
            col = "black", pch = 20)
        par(lty = 1)
        abline(h = 0, v = 0)
        mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$kup, 2]^2 + 
            x$axys[x$kup, 3]^2)
        xray <- c(0, x$axys[x$kup, 2]) * mfac
        yray <- c(0, x$axys[x$kup, 3]) * mfac
        lines(xray, yray)
        mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$jlo, 2]^2 + 
            x$axys[x$jlo, 3]^2)
        xray <- c(0, x$axys[x$jlo, 2]) * mfac
        yray <- c(0, x$axys[x$jlo, 3]) * mfac
        lines(xray, yray)
        par(lty = 3)
        abline(c(0, x$lambda))
        }
    par(lty = 1)
    opar <- par(ask = dev.interactive(orNone = TRUE))
    hist(x$pref[x$axys[, 4] == 1], main = "", xlab = "Preference Score")
    title(main = "Economic Preference Distribution within ICE Wedge", 
        font.main = 3)
    par(opar)
}
`plot.ICEepmap` <-
function (x, xygrid = FALSE, ...) 
{
    if (missing(x)) 
        stop("The first argument to plot.ICEepmap must be an ICEepmap object.")
    if (xygrid == FALSE) {
        xp <- seq(-10, +10, length = 201)
        yp <- xp
        xygrid <- expand.grid(x = xp, y = yp)
    }
    r <- as.vector(sqrt(xygrid$x^2 + xygrid$y^2/x$lambda^2))
    abcos <- abs(xygrid$x - xygrid$y/x$lambda)/r
    xygrid$z <- sign(xygrid$x - xygrid$y/x$lambda) * r^x$beta * 
        abcos^x$gamma
    eta <- x$gamma/x$beta
    contourplot(z ~ x * y, xygrid, cuts = 50, xlab = "Delta Effectiveness (Cost Units)", 
        ylab = "Delta Cost", main = "ICE Economic Preference Map", 
        sub = paste("lambda =", round(x$lambda, digits = 3), 
            ", beta =", round(x$beta, digits = 3), ", gamma=", 
            round(x$gamma, digits = 3), ", eta=", round(eta, 
                digits = 3)), contour = TRUE, labels = FALSE, 
        region = TRUE, col.regions = rainbow(100, end = 0.33), 
        bg = "white")
}
`plot.ICEuncrt` <-
function (x, lfact = 1, swu = FALSE, alibi = FALSE, ...) 
{
    if (missing(x) || !inherits(x, "ICEuncrt")) 
        stop("The first argument to plot(ICEuncrt) must be an ICEuncrt object.")
    if (lfact > 0) 
        lambda <- lfact * x$lambda
    else lambda <- x$lambda
    unit <- x$unit
    if (unit != "cost") 
        unit <- "effe"
    if (swu) {
        if (unit == "cost") { 
            unit <- "effe"
            if (x$lambda != 1) {
                x$t1[1] <- x$t1[1] / x$lambda
                x$t[, 1] <- x$t[, 1] / x$lambda
                x$t1[2] <- x$t1[2] / x$lambda
                x$t[, 2] <- x$t[, 2] / x$lambda
                }
            }
        else {
            unit <- "cost"
            if (x$lambda != 1) {
                x$t1[1] <- x$t1[1] * x$lambda
                x$t[, 1] <- x$t[, 1] * x$lambda
                x$t1[2] <- x$t1[2] * x$lambda
                x$t[, 2] <- x$t[, 2] * x$lambda
                }
            }
        }
    if (lfact != 1) {
        x$lambda <- lambda
        if (unit == "cost") {
            x$t1[1] <- x$t1[1] * lfact
            x$t[, 1] <- x$t[, 1] * lfact
        }
        else {
            x$t1[2] <- x$t1[2] / lfact
            x$t[, 2] <- x$t[, 2] / lfact
        }
    }
    emax <- max(abs(max(x$t[, 1])), abs(min(x$t[, 1])))
    cmax <- max(abs(max(x$t[, 2])), abs(min(x$t[, 2])))
    if (alibi == FALSE) {
        plot(x$t[, 1], x$t[, 2], ann = FALSE, type = "p", ylim = c(-cmax, 
            cmax), xlim = c(-emax, emax))
        par(lty = 1)
        abline(v = 0, h = 0)
        par(lty = 2)
        abline(v = x$t1[1], h = x$t1[2])
        par(lty = 3)
        abline(c(0, 1))
        title(main = paste("ICE Alias Uncertainty for Lambda =", 
            lambda), xlab = "Effectiveness Difference", ylab = "Cost Difference", 
            sub = paste("Units =", unit, ": Bootstrap Reps =", x$R))
        }
    else {
        amax <- max(emax, cmax)
        plot(x$t[, 1], x$t[, 2], ann = FALSE, type = "p", ylim = c(-amax, 
            amax), xlim = c(-amax, amax))
        par(lty = 1)
        abline(v = 0, h= 0)
        par(lty = 2)
        abline(v = x$t1[1], h = x$t1[2])
        par(lty = 3)
        abline(c(0, lambda))
        title(main = paste("ICE Alibi Uncertainty for Lambda =", 
            lambda), xlab = "Effectiveness Difference", ylab = "Cost Difference", 
            sub = paste("Units =", unit, ": Bootstrap Reps =", x$R))
        }
    ICEuncol <- list(df = x$df, lambda = lambda, unit = unit, R = x$R,
        trtm = x$trtm, xeffe = x$xeffe, ycost = x$ycost, effcst = x$effcst,
        t1 = x$t1, t = x$t, seed = x$seed)
    class(ICEuncol) <- "ICEuncrt"
    ICEuncol
}
`plot.ICEwedge` <-
function (x, ...) 
{
    if (missing(x) || !inherits(x, "ICEwedge")) 
        stop("The first argument to plot.ICEwedge must be an ICEwedge object.")
    plot(x$axys[, 2], x$axys[, 3], ann = FALSE, type = "p", pch = 20, 
        col = c("black", "cyan")[unclass(x$axys[, 4]) + 1], bg = "white", 
        ylim = c(-x$ymax, x$ymax), xlim = c(-x$xmax, x$xmax))
    par(lty = 1)
    abline(v = 0)
    abline(h = 0)
    mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$kup, 2]^2 + 
        x$axys[x$kup, 3]^2)
    xray <- c(0, x$axys[x$kup, 2]) * mfac
    yray <- c(0, x$axys[x$kup, 3]) * mfac
    lines(xray, yray)
    mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$jlo, 2]^2 + 
        x$axys[x$jlo, 3]^2)
    xray <- c(0, x$axys[x$jlo, 2]) * mfac
    yray <- c(0, x$axys[x$jlo, 3]) * mfac
    lines(xray, yray)
    par(lty = 2)
    mfac <- 10 * max(x$xmax, x$ymax)/sqrt(x$axys[x$center, 2]^2 + 
        x$axys[x$center, 3]^2)
    xray <- c(0, x$axys[x$center, 2]) * mfac
    yray <- c(0, x$axys[x$center, 3]) * mfac
    lines(xray, yray)
    par(lty = 1)
    title(main = paste("Wedge-Shaped ICE Region with Confidence =", 
        100 * x$conf, "%"), xlab = "Effectiveness Difference", 
        ylab = "Cost Difference", sub = paste("Units =", x$unit, 
            "; lambda =", round(x$lambda, digits = 3), "; Angles =", x$ab))
}
`print.ICEalice` <-
function (x, ...) 
{
    cat("\nICEalice: Acceptability Curves from Bootstrap Uncertainty Distribution...\n")
    cat(paste("\nShadow Price of Health, lambda =", x$lambda))
    cat(paste("\nICE Differences in both Cost and Effectiveness expressed in", 
        x$unit, "units.\n\n"))
    print(x$acc)
    cat("\nICE Quadrant Confidence Level Percentages... (SE, NE, SW, NW)\n")
    print(x$qcl) 
    cat("\n")
}
`print.ICEcolor` <-
function (x, ...) 
{
    cat("\nICEcolor: Economic Preference Coloring of ICE Uncertainty Distribution...\n")
    cat(paste("Shadow Price of Health, lambda:", x$lambda, "\n"))
    cat(paste("Black Points are outside ICE Wedge with Confidence =", 
        100 * x$conf, "%\n"))
    cat(paste("Returns-to-Scale Power, beta:", x$beta, "\n"))
    cat(paste("Preference Shape Power, gamma:", round(x$gamma, 
        digits = 3), "\n"))
    eta <- x$gamma/x$beta
    if (eta > 3 + 2 * sqrt(2)) {
        cat("Power Parameter Ratio, Eta:", round(eta, digits = 3), 
            " > Omega = 3+2*sqrt(2)\n")
        cat("Uncertainty Scatter colored with Highly Directional Preferences lacking Monotonicity.\n\n")
    }
    else if (eta < 3 - 2 * sqrt(2)) {
        cat("Power Parameter Ratio, Eta:", round(eta, digits = 3), 
            " < 1/Omega = 3-2*sqrt(2)\n")
        cat("Uncertainty Scatter colored with Roundish Preferences lacking Monotonicity.\n\n")
    }
    else {
        cat("Power Parameter Ratio, Eta:", round(eta, digits = 3), 
            "\n")
        cat("Uncertainty Scatter colored with Monotonic Preferences & Non-Negative Willingness.\n\n")
    }
    summary(x$pref[x$axys[, 4] == 1])
}
`print.ICEepmap` <-
function (x, ...) 
{
    cat("\nICEepmap: Economic Preference Map...\n")
    cat("Shadow Price of Health, Lambda:", x$lambda, "\n")
    cat("Returns-to-Scale Power, Beta:", x$beta, "\n")
    cat("Preference Shape Power, Gamma:", x$gamma, "\n")
    eta <- x$gamma/x$beta
    if (eta > 3 + 2 * sqrt(2)) {
        cat("Power Parameter Ratio, Eta:", eta, " > Omega = 3+2*sqrt(2)\n")
        cat("This is a Highly Directional ICE map that lacks Monotonicity.\n")
    }
    else if (eta < 3 - 2 * sqrt(2)) {
        cat("Power Parameter Ratio, Eta:", eta, " < 1/Omega = 3-2*sqrt(2)\n")
        cat("This is a Roundish ICE map that lacks Monotonicity.\n")
    }
    else {
        cat("Power Parameter Ratio, Eta:", eta, "\n")
        cat("This is an ICE map with Monotonicity and non-negative Willingness.\n")
    }
}
`print.ICEscale` <-
function (x, ...) 
{
    if (missing(x) || !inherits(x, "ICEscale")) 
        stop("The first argument to print.ICEscale must be an ICEscale object.")
    cat("\nIncremental Cost-Effectiveness (ICE) Lambda Scaling Statistics\n")
    cat(paste("\nSpecified Value of Lambda   =", x$lambda))
    cat(paste("\nCost and Effe Differences are both expressed in", 
        x$unit, "units\n"))
    cat(paste("\nEffectiveness variable Name =", x$xeffe))
    cat(paste("\n     Cost     variable Name =", x$ycost))
    cat(paste("\n  Treatment   factor   Name =", x$trtm))
    cat(paste("\nNew treatment level is =", names(table(x$effcst[, 
        1]))[2], "and Standard level is =", names(table(x$effcst[, 
        1]))[1], "\n"))
    cat(paste("\nObserved  Treatment Diff =", round(x$t1[1], 
        digits = 3)))
    cat(paste("\nStd. Error of Trtm Diff  =", round(x$s1[1], 
        digits = 3), "\n"))
    cat(paste("\nObserved Cost Difference =", round(x$t1[2], 
        digits = 3)))
    cat(paste("\nStd. Error of Cost Diff  =", round(x$s1[2], 
        digits = 3), "\n"))
    cat(paste("\nObserved  ICE  Ratio     =", round(x$t1[2]/x$t1[1], 
        digits = 3), "\n"))
    slam <- x$s1[2]/x$s1[1]
    cat(paste("\nStatistical Shadow Price =", round(slam, digits = 3)))
    cat(paste("\nPower-of-Ten Shadow Price=", 10^(as.integer(log10(slam))), 
        "\n\n"))
}
`print.ICEuncrt` <-
function (x, lfact = 1, swu = FALSE, ...) 
{
    if (missing(x) || !inherits(x, "ICEuncrt")) 
        stop("The first argument to print.ICEuncrt() must be an ICEuncrt object.")
    cat("\nIncremental Cost-Effectiveness (ICE) Bivariate Bootstrap Uncertainty\n")
    if (lfact > 0) 
        lambda <- lfact * x$lambda
    else lambda <- x$lambda
    unit <- x$unit
    if (unit != "cost") 
        unit <- "effe"
    if (swu) {
        if (unit == "cost") { 
            unit <- "effe"
            if (x$lambda != 1) {
                x$t1[1] <- x$t1[1] / x$lambda
                x$t[, 1] <- x$t[, 1] / x$lambda
                x$t1[2] <- x$t1[2] / x$lambda
                x$t[, 2] <- x$t[, 2] / x$lambda
                }
            }
        else {
            unit <- "cost"
            if (x$lambda != 1) {
                x$t1[1] <- x$t1[1] * x$lambda
                x$t[, 1] <- x$t[, 1] * x$lambda
                x$t1[2] <- x$t1[2] * x$lambda
                x$t[, 2] <- x$t[, 2] * x$lambda
                }
            }
        }
    cat(paste("\nShadow Price = Lambda =", lambda))
    cat(paste("\nBootstrap Replications, R =", x$R))
    cat(paste("\nEffectiveness variable Name =", x$xeffe))
    cat(paste("\n     Cost     variable Name =", x$ycost))
    cat(paste("\n  Treatment   factor   Name =", x$trtm))
    cat(paste("\nNew treatment level is =", names(table(x$effcst[, 
        1]))[2], "and Standard level is =", names(table(x$effcst[, 
        1]))[1], "\n"))
    cat(paste("\nCost and Effe Differences are both expressed in", 
        unit, "units\n"))
    if (lfact != 1) {
        x$lambda <- lambda
        if (unit == "cost") {
            x$t1[1] <- x$t1[1] * lfact
            x$t[, 1] <- x$t[, 1] * lfact
        }
        else {
            x$t1[2] <- x$t1[2] / lfact
            x$t[, 2] <- x$t[, 2] / lfact
        }
    }
    cat(paste("\nObserved  Treatment Diff =", round(x$t1[1], 
        digits = 3)))
    cat(paste("\nMean Bootstrap Trtm Diff =", round(mean(x$t[, 
        1]), digits = 3), "\n"))
    cat(paste("\nObserved Cost Difference =", round(x$t1[2], 
        digits = 3)))
    cat(paste("\nMean Bootstrap Cost Diff =", round(mean(x$t[, 
        2]), digits = 3), "\n\n"))
    ICEuncol <- list(df = x$df, lambda = lambda, unit = unit, R = x$R,
        trtm = x$trtm, xeffe = x$xeffe, ycost = x$ycost, effcst = x$effcst,
        t1 = x$t1, t = x$t, seed = x$seed)
    class(ICEuncol) <- "ICEuncrt"
    ICEuncol
}
`print.ICEwedge` <-
function (x, ...) 
{
    cat("\nICEwedge: Incremental Cost-Effectiveness Bootstrap Confidence Wedge...\n")
    cat(paste("\nShadow Price of Health, lambda =", x$lambda))
    cat(paste("\nShadow Price of Health Multiplier, lfact =", 
        x$lfact))
    cat(paste("\nICE Differences in both Cost and Effectiveness expressed in", 
        x$unit, "units."))
    cat(paste("\nICE Angle of the Observed Outcome =", round(x$ia1, 
        digits = 3)))
    cat(paste("\nICE Ratio of the Observed Outcome =", round(x$t1[2]/x$t1[1], 
        digits = 5)))
    cat(paste("\nCount-Outwards  Central ICE Angle Order Statistic =", 
        x$center, "of", x$R))
    cat(paste("\nCounter-Clockwise Upper ICE Angle Order Statistic =", 
        x$kup))
    cat(paste("\nCounter-Clockwise Upper ICE Angle =", round(x$axys[x$kup, 
        1], digits = 3)))
    cat(paste("\nCounter-Clockwise Upper ICE Ratio =", round(x$axys[x$kup, 
        3]/x$axys[x$kup, 2], digits = 5)))
    cat(paste("\n    Clockwise     Lower ICE Angle Order Statistic =", 
        x$jlo))
    cat(paste("\n    Clockwise     Lower ICE Angle =", round(x$axys[x$jlo, 
        1], digits = 3)))
    cat(paste("\n    Clockwise     Lower ICE Ratio =", round(x$axys[x$jlo, 
        3]/x$axys[x$jlo, 2], digits = 5)))
    cat(paste("\nICE Angle Computation Perspective =", x$ab))
    cat(paste("\nConfidence Wedge Subtended ICE Polar Angle =", 
        round(x$subangle, digits = 3)))        
    cat("\n\n")
}