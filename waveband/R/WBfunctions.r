"fill.mat" <- 
function(mat.row, coeffs)
{
    return(coeffs[mat.row[1]:mat.row[2]])
}

"plot.wb" <- 
function(x, col = FALSE, ...)
{
    #
    # Data should be the output from running wave.band
    #
    # Plots the BayesThresh estimate with 99% credible intervals
    # and the data.
    #
    #
    # If col=TRUE, then the plot should look the same as that generated 
    # by wave.band; otherwise a different plot that looks better
    # in black and white is produced.
    #
    n <- length(x$data)
    xtmp <- (1:n)/n
    if(x$param$type != "data")
        y <- test.data(type = x$param$type, rsnr = x$
            param$rsnr, n = n, signal = 1)$y
    if(col) {
        plot(xtmp, x$data, type = "l", xlab = "x", ylab = "y",
            ylim = range(x$data, x$bands$l99, x$
            bands$u99))
        lines(xtmp, x$data, col = 3)
        lines(xtmp, x$bands$l99, col = 2)
        lines(xtmp, x$bands$u99, col = 2)
        lines(xtmp, x$bands$pointest, col = 4)
        if(x$param$type != "data")
            lines(x, y, lty = 2)
    }
    else {
        plot(xtmp, x$data, xlab = "x", ylab = "y", ylim = range(
            x$data, x$bands$l99, x$bands$u99),
            pch = ".")
        lines(xtmp, x$bands$l99)
        lines(xtmp, x$bands$u99)
        lines(xtmp, x$bands$pointest, lwd = 2)
        if(x$param$type != "data")
            lines(x, y, lty = 2)
    }
}

"power.sum" <- 
function(alphas.wd, pow = 2, verbose = TRUE, type = "approx", plotfn = FALSE)
{
    #
    # This function evaluates the sum 
    #
    #   sum_{j,k} alpha_{j,k} psi_{j,k}^{pow}(x)
    #
    # where the psi_{j,k} are mother wavelets of some kind.  The .wd 
    # object alphas.wd contains the alpha_{j,k} in its $D component, 
    # and specifies the wavlet to use in its $filter component.  
    #
    #
    # NB: if pow=2, there should be a coefficient of
    # phi_{0,0}^2(x) to take care of as well.
    #
    # Exact and approximate solutions are available; the approximation
    # is very good and takes about 1/3 the time of the exact solution.
    #
    # Firstly, extract some components of alphas.wd for later use:
    #
    filter.number = alphas.wd$filter$filter.number
    family = alphas.wd$filter$family
    J <- nlevelsWT(alphas.wd)
    #
    # J is log_2(length(data)). Note that this function implicitly 
    # requires J to be at least 5 (data of length 32).
    #
    # Create a .wd object full of zeros for use as a work space, 
    # and a copy that we'll build our estimate in.
    #
    zero.wd <- wd(rep(0, 2^J), filter.number = filter.number, 
        family = family)
    estimate.wd <- zero.wd
    #
    # If asked, do the exact version
    #
    if(type != "approx") {
        if(verbose)
            cat("Evaluating the exact solution\n")
        exact.sum <- rep(0, 2^J)
        if(pow == 2) {
            if(verbose)
                cat("Including overall scaling function\n"
                    )
            phisq <- wr(putC(zero.wd, level = 0, v = 1))^
                pow
            exact.sum <- exact.sum + accessC(alphas.wd,
                level = 0) * phisq
        }
        for(j in 0:(J - 1)) {
            if(verbose)
                cat(c("Starting work on level", j,
                    "\n"))
            for(k in 0:(2^j - 1)) {
                before <- k
                after <- 2^j - k - 1
                psisq <- wr(putD(zero.wd, level = j,
                    v = c(rep(0, before), 1, rep(
                    0, after))))^pow
                exact.sum <- exact.sum + accessD(
                    alphas.wd, level = j)[k + 1] *
                    psisq
            }
        }
        if(type == "exact" && plotfn == TRUE)
            plot(exact.sum, xlab = "x", ylab = "y", type = 
                "l")
        if(type == "exact")
            return(exact.sum)
    }
    #
    # Do levels j = 0,...,J-4 by getting psi^pow_{j,0} and shifting 
    # (if necessary) to get psi^pow_{j,k}.  The approximation is
    # by scaling coefficients at level j+3.
    #
    if(verbose) cat("Evaluating the approximate sum\n")
    if(pow == 2) {
        if(verbose)
            cat("Including overall scaling function\n")
        pow.coeffs <- accessC(wd(wr(putC(zero.wd, level = 0,
            v = 1))^pow, filter.number = filter.number,
            family = family), level = 3)
        tmp <- pow.coeffs * accessC(alphas.wd, level = 0)
        estimate.wd <- putC(estimate.wd, level = 3, v = tmp)
    }
    for(j in 0:(J - 4)) {
        if(verbose)
            cat(c("Starting work on level", j, "\n"))
        pow.coeffs <- accessC(wd(wr(putD(zero.wd, level = j,
            v = c(1, rep(0, 2^j - 1))))^pow, filter.number
             = filter.number, family = family), level = j +
            3)
        pow.coeffs <- rep(pow.coeffs, 2)
        starts <- ((2^j):1) * 8 + 1
        stops <- starts + 2^(j + 3) - 1
        tmp <- t(apply(cbind(starts, stops), 1, fill.mat, 
            pow.coeffs))
        tmp <- accessD(alphas.wd, level = j) * tmp
        tmp <- apply(tmp, 2, sum) + accessC(estimate.wd, level
             = j + 3)
        estimate.wd <- putC(estimate.wd, level = j + 3, v = tmp
            )
    }
    #
    # At finer levels, can only go up as far as level J-1.
    # Use the wavelet coeffs as well as the scaling coeffs to compensate.
    #
    for(j in (J - 3):(J - 1)) {
        if(verbose)
            cat(c("Starting work on level", j, "\n"))
        starts <- ((2^j):1) * (2^(J - 1 - j)) + 1
        stops <- starts + 2^(J - 1) - 1
        psijtothepow.wd <- wd(wr(putD(zero.wd, level = j, v = c(
            1, rep(0, 2^j - 1))))^pow, filter.number = 
            filter.number, family = family)
        #
        # Scaling coeffs first
        #
        pow.coeffs <- rep(accessC(psijtothepow.wd, level = J -
            1), 2)
        tmp <- t(apply(cbind(starts, stops), 1, fill.mat, 
            pow.coeffs))
        tmp <- accessD(alphas.wd, level = j) * tmp
        tmp <- apply(tmp, 2, sum) + accessC(estimate.wd, level
             = J - 1)
        estimate.wd <- putC(estimate.wd, level = J - 1, v = tmp
            )
        #
        # Now the wavelet coefficients
        #
        pow.coeffs <- rep(accessD(psijtothepow.wd, level = J -
            1), 2)
        tmp <- t(apply(cbind(starts, stops), 1, fill.mat, 
            pow.coeffs))
        tmp <- accessD(alphas.wd, level = j) * tmp
        tmp <- apply(tmp, 2, sum) + accessD(estimate.wd, level
             = J - 1)
        estimate.wd <- putD(estimate.wd, level = J - 1, v = tmp
            )
    }
    estimate <- wrow(estimate.wd)
    if(plotfn) {
        if(type == "approx")
            plot(estimate, type = "l", xlab = "x", ylab = 
                "y")
        else {
            plot(exact.sum, type = "l", xlab = "x", ylab = 
                "y", ylim = range(c(estimate, 
                exact.sum)))
            lines(estimate, col = 2)
        }
    }
    if(type == "both")
        return(list(estimate = estimate, exact.sum = exact.sum)
            )
    return(estimate)
}

"test.data" <- 
function(type = "ppoly", n = 512, signal = 1, rsnr = 7, plotfn = FALSE)
{
    x <- seq(0., 1., length = n + 1)[1:n]
    #    elseif strcmp(Name,'Bumps'), 
    #    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    #    hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
    #    wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
    #    sig = zeros(size(t));
    #    for j =1:length(pos)
    #       sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
    #    end 
    if(type == "ppoly") {
        y <- rep(0., n)
        xsv <- (x <= 0.5)
        y[xsv] <- -16. * x[xsv]^3. + 12. * x[xsv]^2.
        xsv <- (x > 0.5) & (x <= 0.75)
        y[xsv] <- (x[xsv] * (16. * x[xsv]^2. - 40. * x[xsv] +
            28.))/3. - 1.5
        xsv <- x > 0.75
        y[xsv] <- (x[xsv] * (16. * x[xsv]^2. - 32. * x[xsv] +
            16.))/3.
    }
    else if(type == "blocks") {
        t <- c(0.10000000000000001, 0.13, 0.14999999999999999,
            0.23000000000000001, 0.25, 0.40000000000000002,
            0.44, 0.65000000000000002, 0.76000000000000001,
            0.78000000000000003, 0.81000000000000005)
        h <- c(4., -5., 3., -4., 5., -4.2000000000000002, 
            2.1000000000000001, 4.2999999999999998, 
            -3.1000000000000001, 2.1000000000000001, 
            -4.2000000000000002
            )
        y <- rep(0., n)
        for(i in seq(1., length(h))) {
            y <- y + (h[i] * (1. + sign(x - t[i])))/2.
        }
    }
    else if(type == "bumps") {
        t <- c(0.10000000000000001, 0.13, 0.14999999999999999,
            0.23000000000000001, 0.25, 0.40000000000000002,
            0.44, 0.65000000000000002, 0.76000000000000001,
            0.78000000000000003, 0.81000000000000005)
        h <- c(4., 5., 3., 4., 5., 4.2000000000000002, 
            2.1000000000000001, 4.2999999999999998, 
            3.1000000000000001, 5.0999999999999996, 
            4.2000000000000002)
        w <- c(0.0050000000000000001, 0.0050000000000000001,
            0.0060000000000000001, 0.01, 0.01, 
            0.029999999999999999, 0.01, 0.01, 
            0.0050000000000000001, 0.0080000000000000002,
            0.0050000000000000001)
        y <- rep(0, n)
        for(j in 1:length(t)) {
            y <- y + h[j]/(1. + abs((x - t[j])/w[j]))^
                4.
        }
    }
    else if(type == "heavi")
        y <- 4. * sin(4. * pi * x) - sign(x - 
            0.29999999999999999) - sign(0.71999999999999997 -
            x)
    else if(type == "doppler") {
        eps <- 0.050000000000000003
        y <- sqrt(x * (1. - x)) * sin((2. * pi * (1. + eps))/
            (x + eps))
    }
    else {
        cat(c("test.data: unknown test function type", type,
            "\n"))
        cat(c("Terminating\n"))
        return("NoType")
    }
    y <- y/sqrt(var(y)) * signal
    ynoise <- y + rnorm(n, 0, signal/rsnr)
    if(plotfn == TRUE) {
        if(type == "ppoly")
            mlab = "Piecewise polynomial"
        if(type == "blocks")
            mlab = "Blocks"
        if(type == "bumps")
            mlab = "Bumps"
        if(type == "heavi")
            mlab = "HeaviSine"
        if(type == "doppler")
            mlab = "Doppler"
        plot(x, y, type = "l", lwd = 2, main = mlab, ylim = 
            range(c(y, ynoise)))
        lines(x, ynoise, col = 2)
        lines(x, y)
    }
    return(list(x = x, y = y, ynoise = ynoise, type = type, rsnr = 
        rsnr))
}

"wave.band" <- 
function(data = 0, alpha = 0.5, beta = 1., filter.number = 8, family = 
    "DaubLeAsymm", bc = "periodic", dev = var, j0 = 3., plotfn = TRUE,
    retvalue = TRUE, n = 128, type = "data", rsnr = 3)
{
    #
    # Data should be a vector of noisy data of length 2^J
    #
    # Alternatively, specify type to be "ppoly", "blocks", "bumps", 
    # "doppler", or "heavi" and wave.band will call test.data.
    #
    if(type == "data") {
        n <- length(data)
        ispow <- !is.na(IsPowerOfTwo(n))
        if(ispow && (n < 1025))
            data <- list(x = (1:n)/n, ynoise = data)
        else {
            cat("Warning: ")
            if(n > 1025)
                cat("data vector is too long (over 1024)\n"
                    )
            else cat("length of data is not a power of two\n"
                    )
            return(NULL)
        }
    }
    if(type != "data") {
        data <- test.data(type = type, rsnr = rsnr, n = n)
        if(is.list(data) == FALSE)
            return(NULL)
    }
    #
    # Estimation of hyperparamters C1 and C2 via universal thresholding;
    # Entirely unchanged from BAYES.THR
    #
    ywd <- wd(data$ynoise, filter.number = filter.number, family = 
        family, bc = bc)
    sigma <- sqrt(dev(accessD(ywd, level = (nlevelsWT(ywd) - 1.))))
    uvt <- threshold(ywd, policy = "universal", type = "soft",
        dev = dev, by.level = FALSE, levels = (nlevelsWT(ywd) - 1.),
        return.threshold = TRUE)
    universal <- threshold(ywd, policy = "manual", value = uvt,
        type = "soft", dev = dev, levels = j0:(nlevelsWT(ywd) -
        1.))
    nsignal <- rep(0., nlevelsWT(ywd))
    sum2 <- rep(0., nlevelsWT(ywd))
    for(j in 0.:(nlevelsWT(ywd) - 1.)) {
        coefthr <- accessD(universal, level = j)
        nsignal[j + 1.] <- sum(abs(coefthr) > 0.)
        if(nsignal[j + 1.] > 0.)
            sum2[j + 1.] <- sum(coefthr[abs(coefthr) > 0.]^
                2.)
    }
    C <- seq(1000., 15000., 50.)
    l <- rep(0., length(C))
    lev <- seq(0., nlevelsWT(ywd) - 1.)
    v <- 2.^( - alpha * lev)
    for(i in 1.:length(C)) {
        l[i] <- 0.5 * sum( - nsignal * (logb(sigma^2. + C[
            i] * v) + 2. * logb(pnorm(( - sigma * sqrt(
            2. * logb(2.^nlevelsWT(ywd))))/sqrt(sigma^2. +
            C[i] * v)))) - sum2/2./(sigma^2. + C[i] * v))
    }
    C1 <- C[l == max(l)]
    tau2 <- C1 * v
    p <- 2. * pnorm(( - sigma * sqrt(2. * logb(2.^nlevelsWT(ywd))))/
        sqrt(sigma^2. + tau2))
    if(beta == 1.)
        C2 <- sum(nsignal/p)/nlevelsWT(ywd)
    else C2 <- (1. - 2.^(1. - beta))/(1. - 2.^((1. - beta) * nlevelsWT(ywd))) * sum(nsignal/p)
    pr <- pmin(1., C2 * 2.^( - beta * lev))
    rat <- tau2/(sigma^2. + tau2)
    #
    # Now we work out the cumulants of the wavelet coefficients.
    #
    K1.wd <- wd(rep(0, n), filter.number = filter.number, family = 
        family, bc = bc)
    K2.wd <- K1.wd
    K3.wd <- K1.wd
    K4.wd <- K1.wd
    bayesthresh.wd <- ywd
    #
    # Deal with the cumulants of the scaling coefficient:
    #
    c00 <- accessC(ywd, level = 0)
    K1.wd <- putC(K1.wd, level = 0, v = c00)
    K2.wd <- putC(K2.wd, level = 0, v = sigma^2)
    #
    # Now the cumulants of the wavelet coefficients:
    #
    for(j in 0.:(nlevelsWT(ywd) - 1.)) {
        coef <- accessD(ywd, level = j)
        w <- ((1. - pr[j + 1.])/pr[j + 1.])/(sqrt((sigma^2 *
            rat[j + 1])/tau2[j + 1.])) * exp(( - rat[j +
            1] * coef^2)/(2 * sigma^2))
        z <- 0.5 * (1. + pmin(w, 1.))
        median <- sign(coef) * pmax(0., rat[j + 1.] * abs(
            coef) - sigma * sqrt(rat[j + 1.]) * qnorm(
            z))
        bayesthresh.wd <- putD(bayesthresh.wd, level = j, v = 
            median)
        g <- 1/(1 + w)
        gt <- 1 - g
        #
        # First cumulant
        #
        k1 <- g * coef * rat[j + 1]
        K1.wd <- putD(K1.wd, level = j, v = k1)
        #
        # Second cumulant
        #
        k2 <- (g * rat[j + 1] * (coef^2 * rat[j + 1] * gt +
            sigma^2))
        K2.wd <- putD(K2.wd, level = j, v = k2)
        #
        # Third cumulant:
        #
        k3 <- (g * gt * coef * rat[j + 1]^2 * (coef^2 * rat[
            j + 1] * (1 - 2 * g) + 3 * sigma^2))
        K3.wd <- putD(K3.wd, level = j, v = k3)
        #
        # Fourth cumulant:
        #
        k4 <- (g * gt * rat[j + 1]^2 * (coef^4 * rat[j + 1]^
            2 * (1 - 6 * g * gt) + 6 * coef^2 * rat[j +
            1] * sigma^2 * (1 - 2 * g) + 3 * sigma^4))
        K4.wd <- putD(K4.wd, level = j, v = k4)
    }
    #
    # Got wavelet coefficient cumulants; now want the cumulants of
    # the function estimate.
    #
    bayesthresh.wr <- wr(bayesthresh.wd)
    K1 <- wrow(K1.wd)
    K2 <- power.sum(K2.wd, pow = 2, verbose = FALSE)
    K3 <- power.sum(K3.wd, pow = 3, verbose = FALSE)
    K4 <- power.sum(K4.wd, pow = 4, verbose = FALSE)
    cumulants <- cbind(K1, K2, K3, K4)
    bands <- t(apply(cumulants, 1, wb.johnson.lims))
    #
    # Now prepare output:
    #
    cumulants <- list(one = cumulants[, 1], two = cumulants[, 2],
        three = cumulants[, 3], four = cumulants[, 4])
    Kr.wd <- list(one = K1.wd, two = K2.wd, three = K3.wd, four = 
        K4.wd)
    bands <- list(pointest = bayesthresh.wr, l80 = bands[, 2],
        u80 = bands[, 3], w80 = bands[, 3] - bands[, 2], l90 = 
        bands[, 4], u90 = bands[, 5], w90 = bands[, 5] - bands[
        , 4], l95 = bands[, 6], u95 = bands[, 7], w95 = bands[
        , 7] - bands[, 6], l99 = bands[, 8], u99 = bands[,
        9], w99 = bands[, 9] - bands[, 8])
    param <- list(alpha = alpha, beta = beta, filter.number = 
        filter.number, family = family, type = type, rsnr = 
        rsnr)
    if(retvalue)
        returnable <- list(data = data$ynoise, cumulants = 
            cumulants, Kr.wd = Kr.wd, bands = bands, param
             = param)
    #
    # And, if asked, draw some pictures:
    #
    if(plotfn == TRUE) {
        plot(data$x, data$ynoise, type = "l", xlab = "x", ylab
             = "y", ylim = range(data$ynoise, bands$l99,
            bands$u99))
        lines(data$x, data$ynoise, col = 3)
        if(type != "data")
            lines(data$x, data$y, lty = 2)
        lines(data$x, bands$l99, col = 2)
        lines(data$x, bands$u99, col = 2)
        lines(data$x, bands$pointest, col = 4)
    }
    if(retvalue)
        return(returnable)
}

"wb.johnson.lims" <- 
function(k, verbose = FALSE)
{
    # Returns the upper and lower 100*siglvl/2 points of the
    # Johnson curve with first four cumulants as given in k, using 
    # algorithms as99 and as100.  This is done for siglvl = 0.2, 0.1, 
    # 0.05, and 0.01.
    #
    # First we calculate the input parameters required by the FORTRAN
    # function "jnsn()":
    #       k[1]: the mean
    #       sd: the standard deviation
    #       rbeta1: sqrt(beta1) in terms of Pearson curves.
    #       NB: it takes the sign of k[3]; this is required
    #         by the algorithm that finds the percentiles.
    #       beta2: parameter beta2 of a Pearson curve;
    #
    sd <- sqrt(k[2])
    rbeta1 <- k[3]/sd^3
    beta2 <- k[4]/sd^4 + 3
    #
    # NB: the "+3" in the anove line is because the usual 
    # definition of beta2 is in terms of moments about the
    # mean; this is what happens when you use cumulants instead
    #
    lims <- rep(0, 9)
    zvalues <- matrix(rep(0, 16), ncol = 2)
    zvalues[, 2] <- c(0.10000000000000001, 0.90000000000000002,
        0.050000000000000003, 0.94999999999999996, 
        0.025000000000000001, 0.97499999999999998, 
        0.0050000000000000001, 0.995)
    zvalues[, 1] <- qnorm(zvalues[, 2])
    #
    # Next use jnsn() to find the Johnson curve with these parameters.
    #
    temp <- .Fortran("jnsn",
        XBAR = as.double(k[1]),
        SD = as.double(sd),
        RB1 = as.double(rbeta1),
        BB2 = as.double(beta2),
        ITYPE = as.integer(0),
        GAMMA = as.double(0),
        DELTA = as.double(0),
        XLAM = as.double(0),
        XI = as.double(0),
        IFAULT = as.integer(0), PACKAGE = "waveband")
    lims[1] <- temp$IFAULT
    #
    # Now use ajv() to get the required points of this distribution, 
    # via the function wb.johnson.lookup.
    #
    lims[2:9] <- apply(zvalues, 1, wb.johnson.lookup, parameters = 
        temp)
    return(lims)
}

"wb.johnson.lookup" <- 
function(zval, parameters)
{
    zvalue <- zval[1]
    temp <- .Fortran("ajv",
        SNV = as.double(zvalue),
        JVAL = as.double(0),
        ITYPE = as.integer(parameters$ITYPE),
        GAMMA = as.double(parameters$GAMMA),
        DELTA = as.double(parameters$DELTA),
        XLAM = as.double(parameters$XLAM),
        XI = as.double(parameters$XI),
        IFAULT = as.integer(0), PACKAGE = "waveband")
    return(temp$JVAL)
}

"wrow" <- 
function(wd, start.level = 0., verbose = FALSE, bc = wd$bc, return.object
     = FALSE, filter.number = wd$filter$filter.number, family = wd$filter$family)
{

    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if(verbose == TRUE)
        cat("Argument checking...")
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- oldClass(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(start.level < 0.)
        stop("start.level must be nonnegative")
    if(start.level >= nlevelsWT(wd))
        stop("start.level must be less than the number of levels"
            )
    if(is.null(wd$filter$filter.number))
        stop("NULL filter.number for wd")
    if(bc != wd$bc)
        warning("Boundary handling is different to original")
    if(wd$type == "station")
        stop("Use convert to generate wst object and then AvBasis or InvBasis"
            )
    if(wd$bc == "interval") {
        warning("All optional arguments ignored for \"wavelets on the interval\" transform"
            )
        return(wr.int(wd))
    }
    type <- wd$type
    filter <- filter.select(filter.number = filter.number, family
         = family)
    LengthH <- length(filter$H)
    if(verbose == TRUE)
        cat("...done\nFirst/last database...")
    r.first.last.c <- wd$fl.dbase$first.last.c[(start.level + 1.):
        (nlevelsWT(wd) + 1.),  ]
    r.first.last.d <- matrix(wd$fl.dbase$first.last.d[(start.level +
        1.):(nlevelsWT(wd)),  ], ncol = 3.)
    ntotal <- r.first.last.c[1., 3.] + r.first.last.c[1., 2.] -
        r.first.last.c[1., 1.] + 1.
    names(ntotal) <- NULL
    C <- wd$C
    #C <- c(rep(0., length = (ntotal - length(C))), C)
    Nlevels <- nlevelsWT(wd) - start.level
    error <- 0.
    if(verbose == TRUE)
        cat("...built\n")
    if(verbose == TRUE) {
        cat("Reconstruction...")
        error <- 1.
    }
    ntype <- switch(type,
        wavelet = 1.,
        station = 2.)
    if(is.null(ntype))
        stop("Unknown type of decomposition")
    nbc <- switch(bc,
        periodic = 1.,
        symmetric = 2.)
    if(is.null(nbc))
        stop("Unknown boundary handling")
    if(!is.complex(wd$D)) {
        wavelet.reconstruction <- .C("wavereconsow",
            C = as.double(C),
            D = as.double(wd$D),
            H = as.double(filter$H),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(Nlevels),
            firstC = as.integer(r.first.last.c[, 1.]),
                lastC = as.integer(r.first.last.c[, 2.]),
            offsetC = as.integer(r.first.last.c[, 3.]),
                firstD = as.integer(r.first.last.d[, 1.]),
            lastD = as.integer(r.first.last.d[, 2.]),
            offsetD = as.integer(r.first.last.d[, 3.]),
                type = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE="waveband")
    }
    else {
        wavelet.reconstruction <- .C("comwr",
            CR = as.double(Re(C)),
            CI = as.double(Im(C)),
            LengthC = as.integer(length(C)),
            DR = as.double(Re(wd$D)),
            DI = as.double(Im(wd$D)),
            LengthD = as.integer(length(wd$D)),
            HR = as.double(Re(filter$H)),
            HI = as.double(Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double(Im(filter$G)),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(Nlevels),
            firstC = as.integer(r.first.last.c[, 1.]),

                lastC = as.integer(r.first.last.c[
                , 2.]),
            offsetC = as.integer(r.first.last.c[, 3.]),

                firstD = as.integer(r.first.last.d[
                , 1.]),
            lastD = as.integer(r.first.last.d[, 2.]),
            offsetD = as.integer(r.first.last.d[, 3.]),

                ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE="waveband")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.reconstruction$error
    if(error != 0.) {
        cat("Error code returned from wavereconsow: ", error,
            "\n")
        stop("wavereconsow returned error")
    }
    fl.dbase <- wd$fl.dbase
    if(!is.complex(wd$D)) {
        l <- list(C = wavelet.reconstruction$C, D = 
            wavelet.reconstruction$D, fl.dbase = fl.dbase,
            nlevels = nlevelsWT(wd), filter = filter, type = 
            type, bc = bc, date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.reconstruction$
            CR, imaginary = wavelet.reconstruction$CI), D = 
            complex(real = wavelet.reconstruction$DR, imaginary = 
            wavelet.reconstruction$DI), fl.dbase = fl.dbase,
            nlevels = nlevelsWT(wd), filter = filter, type = 
            type, bc = bc, date = date())
    }
    oldClass(l) <- "wd"
    if(return.object == TRUE)
        return(l)
    else return(accessC(l))
    stop("Shouldn't get here\n")
}


