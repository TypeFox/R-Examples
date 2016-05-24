squaredgain.wt.filter <- function(filter, level = 1, N = NULL, draw.bands = TRUE, wavelet = TRUE) {
    if(class(filter) == "numeric") {
        filter <- wt.filter(filter)
    }
    if(class(filter) == "character") {
        filter <- wt.filter(tolower(filter), level = level)
    }
    if(class(filter) != "wt.filter") {
        stop("Invalid argument: filter must be of class vector, character, or string.")
    }
    if(wavelet) {
        filter.vector <- filter@h
        wavelet.string <- "Wavelet Filter"
    }
    else {
        filter.vector <- filter@g
        wavelet.string <- "Scaling Filter"
    }

    if(length(filter.vector) < 1024 && is.null(N)) {
        N <- 1024
    }
    if(length(filter.vector) > 1024 && is.null(N)) {
        j <- ceiling(log(length(filter.vector))/log(2))
        N <- 2^j
    }
 
    j <- filter@level

    if(N-length(filter.vector) > 0) {
        zeroes <- rep(0, N-length(filter.vector))
        filter.vector <- c(filter.vector, zeroes)
    }

    Hf <- fft(filter.vector)
    absHf <- abs(Hf)
    squaredHf <- absHf^2
    x <- seq((1/N),1/2, (1/N))

    if(filter@wt.name != "haar") {
        if(filter@wt.name == "none") {
            filter.name <- ""
        }
        else {
            filter.name <- toupper(filter@wt.name)
        }
    }
    else {
        filter.name <- "Haar"
    }

    plot(x, squaredHf[1:(N*.5)], type = "l", xlab = "f", ylab = "", main = paste("Squared gain function for ", filter.name, "Level",j,wavelet.string))

    if(draw.bands) {
        abline(v = 1/(2^(j+1)), lty = 2)
        abline(v = 1/(2^j), lty = 2)
    }
}

