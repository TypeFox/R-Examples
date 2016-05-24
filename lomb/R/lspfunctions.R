lsp <- function(x, times = NULL, from = NULL, to = NULL, type = c("frequency", "period"), ofac = 1, alpha = 0.01, plot = TRUE, ...) {    
    type <- match.arg(type)
    
    if (ofac != floor(ofac)) {
        ofac <- floor(ofac)
        warning("ofac coerced to integer")
    }
    if (ofac < 1) {
        ofac <- 1
        warning("ofac must be integer >=1. Set to 1")
    }
    
    
    if (!is.null(times)) {
        if (!is.vector(times)) 
            stop("no multivariate methods available")
        if (length(x) != length(times)) 
            stop("Length of data and times vector must be equal")
        names <- c(deparse(substitute(times)), deparse(substitute(x)))
    }
    
    if (is.null(times) && is.null(ncol(x))) {
        names <- c("Time", deparse(substitute(x)))
        times <- 1:length(x)
    }
    
    if (is.matrix(x) || is.data.frame(x)) {
        
        if (ncol(x) > 2) 
            stop("no multivariate methods available")
        
        if (ncol(x) == 2) {
            names <- colnames(x)
            times <- x[, 1]
            x <- x[, 2]
        }
    }
    
    times <- times[!is.na(x)]
    x <- x[!is.na(x)]
    
    nobs <- length(x)
    if (nobs < 2) 
        stop("time series must have at least two observations")
    
    times <- as.numeric(times)
    start <- min(times)
    end <- max(times)
    av.int <- mean(diff(times))
    
    o <- order(times)
    times <- times[o]
    x <- x[o]
    
    y <- cbind(times, x)
    colnames(y) <- names
    
    datanames <- colnames(y)
    t <- y[, 1]
    y <- y[, 2]
    
    n <- length(y)
    tspan <- t[n] - t[1]
    fr.d <- 1/tspan
    step <- 1/(tspan * ofac)
    
    if (type == "period") {
        hold <- from
        from <- to
        to <- hold
        if (!is.null(from)) 
            from <- 1/from
        if (!is.null(to)) 
            to <- 1/to
    }
    
    if (is.null(to)) {
        f.max <- floor(0.5 * n * ofac) * step
    } else {
        f.max <- to
    }
    
    freq <- seq(fr.d, f.max, by = step)
    if (!is.null(from)) 
        freq <- freq[freq >= from]
    n.out <- length(freq)
    if (n.out == 0) 
        stop("erroneous frequency range specified ")
    
    x <- t * 2 * pi
    y <- y - mean(y)
    
    norm <- 1/(2 * var(y))
    
    w <- 2 * pi * freq
    PN <- rep(0, n.out)
    for (i in 1:n.out) {
        wi <- w[i]
        tau <- 0.5 * atan2(sum(sin(wi * t)), sum(cos(wi * t)))/wi
        arg <- wi * (t - tau)
        cs <- cos(arg)
        sn <- sin(arg)
        A <- (sum(y * cs))^2
        B <- sum(cs * cs)
        C <- (sum(y * sn))^2
        D <- sum(sn * sn)
        PN[i] <- A/B + C/D
    }
    
    PN <- norm * PN
    
    PN.max <- max(PN)
    peak.freq <- freq[PN == PN.max]
    if (type == "period") 
        peak.at <- c(1/peak.freq, peak.freq) else peak.at <- c(peak.freq, 1/peak.freq)
    
    scanned <- if (type == "frequency") 
        freq else 1/freq
    if (type == "period") {
        scanned <- scanned[n.out:1]
        PN <- PN[n.out:1]
    }
    effm <- 2 * n.out/ofac
    level <- -log(1 - (1 - alpha)^(1/effm))
    exPN <- exp(-PN.max)
    p <- effm * exPN
    if (p > 0.01) 
        p <- 1 - (1 - exPN)^effm
    
    sp.out <- list(scanned = scanned, power = PN, data = datanames, n = n, type = type, ofac = ofac, n.out = n.out, alpha = alpha, sig.level = level, 
        peak = PN.max, peak.at = peak.at, p.value = p)
    
    class(sp.out) <- "lsp"
    
    if (plot) {
        plot(sp.out, ...)
        return(invisible(sp.out))
    } else return(sp.out)
}



plot.lsp <- function(x, type = "l", main = "Lomb-Scargle Periodogram", xlab = NULL, ylab = "normalised power", level = TRUE, log = NULL, ...) {    

    if (is.null(xlab)) 
        xlab <- x$type
    if (is.null(log)) {
        if (x$type == "period") 
            log <- "x" else log <- ""
    }
    
    ymin <- 0
    if (!is.null(x$sig.level)) {
        ymax <- max(x$sig.level, max(x$power)) * 1.25
    } else {
        ymax <- max(x$power) * 1.25
    }
    
    plot(range(x$scanned), c(ymin, ymax), log = log, type = "n", main = main, xlab = xlab, ylab = ylab, ...)
    lines(x$scanned, x$power, type = type)
    
    if (level == TRUE) {
        if (!is.null(x$sig.level)) 
            abline(h = x$sig.level, lty = 2)
    }
    return(invisible(x))
}


summary.lsp <- function(object,...) {
    first <- object$type
    if (first == "frequency") {
        second <- "At period"
    } else {
        second <- "At frequency"
    }
    first <- paste("At ", first)
    from <- min(object$scanned)
    to <- max(object$scanned)
    
    Value <- c(object$data[[1]], object$data[[2]], object$n, object$type, object$ofac, from, to, object$n.out, object$peak, object$peak.at[[1]], object$peak.at[[2]], object$p.value)
    options(warn = -1)
    for (i in 1:length(Value)) {
        if (!is.na(as.numeric(Value[i]))) 
            Value[i] <- format(as.numeric(Value[i]), digits = 5)
    }
    options(warn = 0)
    nmes <- c("Time", "Data", "n", "Type", "Oversampling", "From", "To", "# frequencies", "PNmax", first, second, "P-value (PNmax)")
    report <- data.frame(Value, row.names = nmes)
    report
}


randlsp <- function(repeats=1000, x, times = NULL, from = NULL, to = NULL, type = c("frequency", "period"), ofac = 1, alpha = 0.01, plot = TRUE, trace = TRUE, ...) {
	
	if (is.ts(x)){
		x=as.vector(x)
	}
    
    if (!is.vector(x)) {
        times <- x[, 1]
        x <- x[, 2]
    }
    
    if (plot==TRUE){
    	op <- par(mfrow = c(2, 1))
    }
    
    realres <- lsp(x, times, from, to, type, ofac, alpha, plot = plot, ...)
    realpeak <- realres$peak
    pks <- NULL
    if (trace == TRUE) 
        cat("Repeats: ")
    for (i in 1:repeats) {
        randx <- sample(x, length(x))  # scramble data sequence
        randres <- lsp(randx, times, from, to, type, ofac, alpha, plot = F)
        pks <- c(pks, randres$peak)
        if (trace == TRUE) {
            if (i/10 == floor(i/10)) 
                cat(i, " ")
        }
    }
    if (trace == TRUE) 
        cat("\n")
    
    prop <- length(which(pks >= realpeak))
    p.value <- prop/repeats
    
    if (plot == TRUE) {
    	mx=max(c(pks,realpeak))*1.25
        hist(pks, xlab = "Peak Amplitude", xlim=c(0,mx), main=paste("P-value: ",p.value))
        abline(v = realpeak)
        par(op)
    }
    
    res=realres[-(8:9)]
    res=res[-length(res)]
    res$random.peaks = pks
    res$repeats=repeats
    res$p.value=p.value
    class(res)="randlsp"
    return(invisible(res))
}

summary.randlsp <- function(object,...) {
    first <- object$type
    if (first == "frequency") {
        second <- "At period"
    } else {
        second <- "At frequency"
    }
    first <- paste("At ", first)
    from <- min(object$scanned)
    to <- max(object$scanned)
    
    Value <- c(object$data[[1]], object$data[[2]], object$n, object$type, object$ofac, from, to, object$n.out, object$peak, object$peak.at[[1]], object$peak.at[[2]], object$repeats,object$p.value)
    options(warn = -1)
    for (i in 1:length(Value)) {
        if (!is.na(as.numeric(Value[i]))) 
            Value[i] <- format(as.numeric(Value[i]), digits = 5)
    }
    options(warn = 0)
    nmes <- c("Time", "Data", "n", "Type", "Oversampling", "From", "To", "# frequencies", "PNmax", first, second, "Repeats","P-value (PNmax)")
    report <- data.frame(Value, row.names = nmes)
    report
}

