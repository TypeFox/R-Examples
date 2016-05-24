tar <-
function (y, p1, p2, d, is.constant1 = TRUE, is.constant2 = TRUE, 
    transform = "no", center = FALSE, standard = FALSE, estimate.thd = TRUE, 
    threshold, method=c("MAIC","CLS")[1], a = 0.05, b = 0.95, 
    order.select = TRUE, print = FALSE)
{


cvar <-
function (x, df = 1) 
{
    # x is assumed to be of zero expectation 
    sum(x^2)/(length(x) - df)
}

fmaic <-
function (R, n, order.select = FALSE) 
{
    k <- dim(R)[2]
    v <- R[-1, k]^2
    v <- rev(cumsum(rev(v)))
    k <- k - 1
    AIC <- n * log(v/n) + 2 * seq(k)
    order <- (1:k)[AIC == min(AIC)]
    like <- (-n * log(v/n))/2 - n/2 - n * log(2 * 3.14159)/2
    if (order.select) {
        maic <- min(AIC)
        return(list(MAIC = maic, order = ((1:k)[AIC == maic])[1] - 
            1, AIC = AIC, like = like[order]))
    }
    else return(list(MAIC = AIC[length(AIC)], order = k - 1, 
        AIC = AIC, like = like[length(AIC)]))

}
cat1 <-
function (..., print = TRUE, file = "", sep = " ", fill = FALSE, 
    labels = NULL, append = FALSE) 
{
    if (print) {
        if (is.character(file)) 
            if (file == "") 
                file <- stdout()
            else if (substring(file, 1, 1) == "|") {
                file <- pipe(substring(file, 2), "w")
                on.exit(close(file))
            }
            else {
                file <- file(file, ifelse(append, "a", "w"))
                on.exit(close(file))
            }
         cat(..., file=file, sep=sep, fill=fill, labels=labels, append=append)
    }
    invisible(NULL)
}

makedata <-
function (dataf, p1, p2, d, is.constant1 = TRUE, is.constant2 = TRUE, 
    thd.by.phase = FALSE) 
{
    n <- dim(dataf)[1]
    nseries <- dim(dataf)[2]
    start <- max(c(p1, p2, d)) + 1
    p <- max(c(p1, p2))
    xy <- NULL
    for (i in (1:nseries)) xy <- rbind(xy, cbind(makexy(dataf[, 
        i], p, start, d, thd.by.phase = thd.by.phase), i))
    xy <- dna(xy)
    xy <- xy[s <- order(xy[, p + 2]), ]
    xy1 <- setxy(xy, p, p1, nseries, is.constant1)
    xy2 <- setxy(xy, p, p2, nseries, is.constant2)
    list(xy1 = xy1, xy2 = xy2, sort.list = s)
}

makexy <-
function (x, p, start, d, thd.by.phase = FALSE) 
{
    n <- length(x)
    xy <- NULL
    for (i in (1:p)) xy <- cbind(xy, x[(start - i):(n - i)])
    if (thd.by.phase) 
        xy <- cbind(xy, x[start:n], x[(start - 1):(n - 1)] - 
            x[(start - 2):(n - 2)])
    else xy <- cbind(xy, x[start:n], x[(start - d):(n - d)])
    xy
}

setxy<-
function (old.xy, p, p1, nseries, is.coefficient = TRUE) 
{
    if (p1 >= 1) 
        s <- 1:p1
    else s <- NULL
    n <- dim(old.xy)[1]
    new.xy <- old.xy[, c(s, (p + 1):(p + 3)), drop = FALSE]
    temp <- new.xy[, s, drop = FALSE]
    if (is.coefficient) 
        new.xy <- cbind(1, new.xy)
    if (is.coefficient) 
        temp <- cbind(rep(1, n), temp)
    if (nseries == 1) 
        return(new.xy)
    for (i in rev(2:nseries)) {
        select <- old.xy[, p + 3] == i
        zero <- 0 * temp
        zero[select, ] <- temp[select, ]
        new.xy <- cbind(zero, new.xy)
    }
    new.xy
}

findstart <-
function (x, nseries, indexid, p) 
{
    m <- dim(x)[1]
    amax <- 0
    for (i in (1:nseries)) {
        amax <- max(amax, (1:m)[(cumsum(x[, indexid] == i) == 
            p)])
    }
    amax
}

dna <-
function (x) 
{
    x[!apply(x, 1, any.na), ]
}

any.na <-
function (x) 
{
    any(x == "NA")
}

revm <-
function (m) 
{
    apply(m, 2, rev)
}


if(method=="MAIC") MAIC=TRUE else MAIC=FALSE 

dataf=y

            aic.no.thd <- NA
    if (!is.matrix(dataf)) {
        temp <- cbind(dataf, dataf)
        dataf <- temp[, 1, drop = FALSE]
    }
        dataf <- switch(transform, log = log(dataf), log10 = log10(dataf), 
        sqrt = sqrt(dataf), no = dataf)
    means <- apply(dataf, 2, mean, na.rm = TRUE)
    stds <- apply(dataf, 2, var, na.rm = TRUE)^0.5
    if (standard) 
        dataf <- apply(dataf, 2, standard)
    nseries <- dim(dataf)[2]
    res <- makedata(dataf, p1, p2, d, is.constant1, is.constant2, 
        thd.by.phase = FALSE)
    xy1 <- res$xy1
    xy2 <- res$xy2
    sort.l <- res$sort.list
    m <- dim(xy1)[1]
    q1 <- dim(xy1)[2]
    q2 <- dim(xy2)[2]
    s <- (q1 - 3 - p1):(q1 - 3)
    temp <- xy1[, s]
    xy1 <- cbind(temp, xy1[, -s])
    lab <- deparse(match.call()$y)
    if (p1 >= 1) 
        sublab <- paste("lag", 1:p1, sep = "")
    else sublab <- NULL
    if (is.null(lab)) 
        lab <- ""
    if (lab == "") 
        dimnames(xy1)[[2]] <- as.vector(c(outer(c("intercept", 
            sublab), lab, paste, sep = " "), "", "", ""))
    else dimnames(xy1)[[2]] <- as.vector(c(outer(c("intercept", 
        sublab), lab, paste, sep = "-"), "", "", ""))
    s <- (q2 - 3 - p2):(q2 - 3)
    temp <- xy2[, s]
    xy2 <- cbind(temp, xy2[, -s])
    if (p2 >= 1) 
        sublab <- paste("lag", 1:p2, sep = "")
    else sublab <- NULL
    if (lab == "") 
        dimnames(xy2)[[2]] <- as.vector(c(outer(c("intercept", 
            sublab), lab, paste, sep = " "), "", "", ""))
    else dimnames(xy2)[[2]] <- as.vector(c(outer(c("intercept", 
        sublab), lab, paste, sep = "-"), "", "", ""))
    aic1 <- aic2 <- rss1 <- rss2 <- rep(10^10, m)
    like1 <- like2 <- rep(-10^10, m)
    lbound <- sum(dataf[dataf != "NA"] == min(dataf[dataf != 
        "NA"]))
    ubound <- sum(dataf[dataf != "NA"] == max(dataf[dataf != 
        "NA"]))
    i1 <- max(c(q1 - 3, lbound + 1, 2 * p1 + 1, d, findstart(xy1, 
        nseries, q1, p1 + 2)))
    i1 <- max(i1, floor(a * m))
    i2 <- m - max(c(q2 - 3, ubound + 1, 2 * p2 + 1, d, findstart(revm(xy2), 
        nseries, q2, p2 + 2))) - 1
    i2 <- min(i2, ceiling(b * m))
    s <- -((q1 - 1):q1)
    R <- qr.R(qr(xy1[1:i1, s]))
    if (estimate.thd) {
        posy <- q1 - 2
        rss1[i1] <- (R[posy, posy])^2
        res.fmaic <- fmaic(R, i1, order.select = order.select)
        like1[i1] <- res.fmaic$like
        aic1[i1] <- res.fmaic$MAIC
        for (i in ((i1 + 1):i2)) {
            R <- qr.R(qr(rbind(R, xy1[i, s])))
            rss1[i] <- (R[posy, posy])^2
            res.fmaic <- fmaic(R, i, order.select = order.select)
            like1[i] <- res.fmaic$like
            aic1[i] <- res.fmaic$MAIC
        }
        s <- -((q2 - 1):q2)
        posy <- q2 - 2
        R <- qr.R(qr(xy2[(i2 + 1):m, s]))
        rss2[i2] <- (R[posy, posy])^2
        res.fmaic <- fmaic(R, m - i2, order.select = order.select)
        like2[i2] <- res.fmaic$like
        aic2[i2] <- res.fmaic$MAIC
        for (i in rev(i1:(i2 - 1))) {
            R <- qr.R(qr(rbind(R, xy2[i + 1, s])))
            rss2[i] <- (R[posy, posy])^2
            res.fmaic <- fmaic(R, m - i, order.select = order.select)
            like2[i] <- res.fmaic$like
            aic2[i] <- res.fmaic$MAIC
        }
        rss <- rss1 + rss2
        thdindex <- ((1:m)[rss == min(rss)])[1]
        aic <- aic1 + aic2
        if (MAIC) {
            aic[-(i1:i2)] <- 10^10
            thdindex <- ((1:m)[aic == (aic.no.thd <- min(aic))])[1]
        }
        thd <- xy1[thdindex, q1 - 1]
    }
    else {
        thd <- threshold
        thdindex <- ((1:m)[xy1[, q1 - 1] > thd])[1] - 1
    }
    s <- -((q1 - 1):q1)
    qr1 <- qr(xy1[1:thdindex, s])
    R1 <- qr.R(qr1)
    if(MAIC) {order1=fmaic(R1, thdindex, order.select = order.select)$order+1
subset1=1:order1
    p1=order1-is.constant1} 
else subset1=-((q1 - 2):q1)

        qr1 <- lsfit(x.regime1 <- xy1[1:thdindex, subset1, 
        drop = FALSE], y.regime1 <- xy1[1:thdindex, q1 - 2], intercept = FALSE)
dxy1=xy1[, subset1, drop = FALSE] # for diagnostics
dxy1[-(1:thdindex),]=0
    s <- -((q2 - 1):q2)
    qr2 <- qr(xy2[(thdindex + 1):m, s])
    R2 <- qr.R(qr2)
if(MAIC) {order2=fmaic(R2, m-thdindex, order.select = order.select)$order+1
subset2=1:order2
    p2=order2-is.constant1} 
else subset2=-((q2 - 2):q2)
        qr2 <- lsfit(x.regime2 <- xy2[-(1:thdindex), subset2, 
        drop = FALSE], y.regime2 <- xy2[-(1:thdindex), q2 - 2], intercept = FALSE)
dxy2=xy2[,subset2, drop = FALSE] # for diagnostics
dxy2[(1:thdindex),]=0
    cat1(print = print, "time series included in this analysis is: ", 
        lab, "\n")
    cat1(print = print, "SETAR(2,", p1, ",", p2, ") model")
               cat1(print = print, " delay =", d, "\n")
        if (MAIC) 
        method <- "Minimum AIC"
    else method <- "CLS"
    cat1(print = print, "estimated threshold = ", signif(thd, 
        4), " from a", method, " fit with thresholds \nsearched from the ", 
        round(i1/m * 100), " percentile to the  ", round(i2/m * 
            100), " percentile of all data.\n")
    cat1(print = print, "The estimated threshold is the ", 
        signif(thdindex/m * 100, 3), " percentile of\nall data.\n")    
    cat1(print = print, "lower regime: \n")
    if (print) 
        ls.print(qr1)
    cat1(print = print, "\n\n (unbiased) RMS \n")
    rms1 <- tapply(qr1$residuals, xy1[1:thdindex, q1], cvar, 
        df = p1 + is.constant1)
    cat1(print = print,  signif(rms1, 4), "\n")
    cat1(print = print, " with no of data falling in the regime being \n")
    n1 <- tapply(qr1$residuals, xy1[1:thdindex, q1], length)
    cat1(print = print, rbind(lab, signif(n1, 4)), "\n")
    cat1(print = print, "\n\n (max. likelihood) RMS for each series (denominator=sample size in the regime) \n")
    cat1(print = print, rbind(lab, signif((rms1 * (n1 - 
        p1 - 1))/n1, 4)), "\n")
        cat1(print = print, "\n\n upper regime: \n")
    if (print) 
        ls.print(qr2)
    cat1(print = print, "\n\n (unbiased) RMS \n")
    rms2 <- tapply(qr2$residuals, xy2[(thdindex + 1):m, q2], 
        cvar, df = p2 + is.constant2)
    cat1(print = print,  signif(rms2, 4), "\n")
    cat1(print = print, " with no of data falling in the regime being \n")
    n2 <- tapply(qr2$residuals, xy2[(thdindex + 1):m, q2], length)
    cat1(print = print, signif(n2, 4), "\n")
    cat1(print = print, "\n\n (max. likelihood) RMS for each series (denominator=sample size in the regime)\n")
    cat1(print = print, signif((rms2 * (n2 - 
        p2 - 1))/n2, 4), "\n")
    AIC <- signif(n1 * log((rms1 * (n1 - p1 - is.constant1))/n1) + n2 * 
        log((rms2 * (n2 - p2 - is.constant2))/n2) + 
        +(n1+n2)*(1+log(2*pi))+    2 * (p1 + p2 + is.constant1+is.constant2 + 1), 
        4)
    cat1(print = print, "\n Nominal AIC is ", AIC, "\n")
    residuals <- c(qr1$residuals, qr2$residuals)
    residuals[sort.l] <- residuals
        std.res <- c(qr1$residuals/as.double(rms1^0.5), qr2$residuals/as.double(rms2^0.5))
    std.res[sort.l] <- std.res
dxy1[sort.l,]=dxy1
dxy2[sort.l,]=dxy2

                     
    res=list(dxy1=dxy1,dxy2=dxy2,p1=p1,q1=q1,d=d,qr1=qr1,qr2=qr2,x.regime1 = x.regime1, y.regime1 = y.regime1, 
        x.regime2 = x.regime2, y.regime2 = y.regime2, thd = thd, 
        thdindex = thdindex, qr1 = qr1, qr2 = qr2, i1 = i1, i2 = i2, 
        x = xy1[, q1 - 1], m = m, rss1 = as.vector(rms1 * (n1 - 
            p1 - 1)), rss2 = as.vector(rms2 * (n2 - p2 - 2)), 
        n1 = as.vector(n1), n2 = as.vector(n2), std.res = std.res,
        p1=p1, p2=p2, rms1=rms1,rms2=rms2, 
        is.constant1=is.constant1,is.constant2=is.constant2,
        residuals = residuals, AIC = AIC, aic.no.thd = aic.no.thd,y=y, 
        like = max(like1 + like2),method=method)
class(res)="TAR"
invisible(res)
}
