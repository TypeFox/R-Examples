"spec.ls" <-
function (x, y=NULL, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
    fast = TRUE, type = "lomb", demean = FALSE, detrend = TRUE, 
    plot.it = TRUE, na.action = na.fail, ...) 
{
#    series <- deparse(substitute(x))
    if (NCOL(x)==2){
      series <- deparse(substitute(x))
      ti <- x[,1]
      x <- x[,2]}
#    series <- deparse(substitute(y))}
    else {
       series <- deparse(substitute(y))
      ti <- x
      x <- y
    }
    x <- na.action(as.ts(x))
    xfreq <- frequency(x)
    x <- as.matrix(x)
    N <- nrow(x)
    nser <- ncol(x)
    if (!is.null(spans)) 
        if (is.tskernel(spans)) 
            kernel <- spans
        else kernel <- kernel("modified.daniell", spans%/%2)
    if (!is.null(kernel) && !is.tskernel(kernel)) 
        stop("must specify spans or a valid kernel")
    if (detrend) {
        t <- ti - mean(ti)
        sumt2 <- sum(t^2)
        for (i in 1:ncol(x)) x[, i] <- x[, i] - mean(x[, i]) - 
            sum(x[, i] * t) * t/sumt2
    }
    else if (demean) {
        x <- sweep(x, 2, colMeans(x))
    }
    x <- spec.taper(x, taper)
    u2 <- (1 - (5/8) * taper * 2)
    u4 <- (1 - (93/128) * taper * 2)
    if (pad > 0) {
        x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
        N <- nrow(x)
    }
    NewN <- if (fast) 
        nextn(N)
    else N
    x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
    N <- nrow(x)
    Nspec <- floor(N/2)
    freq <- seq(from = xfreq/N, by = xfreq/N, length = Nspec)
    pgram <- array(NA, dim = c(N, ncol(x), ncol(x)))
    freq.temp <- seq(from = xfreq/N, by = xfreq/N, length = N)
    if (type == "lomb") {
        for (i in 1:ncol(x)) {
            for (j in 1:ncol(x)) {
                for (k in 1:length(freq.temp)) {
                  tao <- atan(sum(sin(2 * pi * freq.temp[k] * 
                                      ti))/sum(cos(2 * pi * freq.temp[k] * ti)))/(2 * freq.temp[k])
                  pgram[k, i, j] <- 0.5 * ((sum(x[1:length(ti)]* cos(2 * pi * freq.temp[k] * (ti - tao))))^2/sum((cos(2 * 
                                                                                                         pi * freq.temp[k] * (ti - tao)))^2) + (sum(x[1:length(ti)] *  sin(2 * pi * freq.temp[k] * (ti - tao))))^2/sum((sin(2 * pi * freq.temp[k] * (ti - tao)))^2))
                  pgram[1, i, j] <- 0.5 * (pgram[2, i, j] + pgram[N, 
                                                                  i, j])
                }
              }
          }
      }
    if (type == "ft") {
        for (i in 1:ncol(x)) {
            for (j in 1:ncol(x)) {
                for (k in 1:length(freq.temp)) pgram[k, i, j] <- ((sum(x * 
                  cos(2 * pi * freq.temp[k] * ti)))^2 + (sum(x * 
                  sin(2 * pi * freq.temp[k] * ti)))^2)/N
            }
        }
    }
    if (!is.null(kernel)) {
        for (i in 1:ncol(x)) for (j in 1:ncol(x)) pgram[, i, 
            j] <- kernapply(pgram[, i, j], kernel, circular = TRUE)
        df <- df.kernel(kernel)/(u4/u2^2)
        bandwidth <- bandwidth.kernel(kernel) * xfreq/N
    }
    else {
        df <- 2/(u4/u2^2)
        bandwidth <- sqrt(1/12) * xfreq/N
    }
    pgram <- pgram[1:Nspec, , , drop = FALSE]
    spec <- matrix(NA, nrow = Nspec, ncol = nser)
    for (i in 1:nser) spec[, i] <- Re(pgram[1:Nspec, i, i])
    if (nser == 1) {
        coh <- phase <- NULL
    }
    else {
        coh <- phase <- matrix(NA, nrow = Nspec, ncol = nser * 
            (nser - 1)/2)
        for (i in 1:(nser - 1)) {
            for (j in (i + 1):nser) {
                coh[, i + (j - 1) * (j - 2)/2] <- Mod(pgram[, 
                  i, j])^2/(spec[, i] * spec[, j])
                phase[, i + (j - 1) * (j - 2)/2] <- Arg(pgram[, 
                  i, j])
            }
        }
    }
    for (i in 1:nser) spec[, i] <- spec[, i]/u2
    spec <- drop(spec)
    spg.out <- list(freq = freq, spec = spec, coh = coh, phase = phase, 
        kernel = kernel, df = df, bandwidth = bandwidth, n.used = nrow(x), 
        series = series, snames = colnames(x), method = ifelse(!is.null(kernel), 
            "Smoothed Lomb-Scargle Periodogram", "Raw Lomb-Scargle Periodogram"), 
        taper = taper, pad = pad, detrend = detrend, demean = demean)
    class(spg.out) <- "spec"
    if (plot.it) {
        plotSpecLs(spg.out, ...)
        return(invisible(spg.out))
    }
    else return(spg.out)
}
