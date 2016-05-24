spectrum0 <- function (x, max.freq = 0.5, order = 1, max.length = 200) 
{
    x <- as.matrix(x)
    if (!is.null(max.length) && nrow(x) > max.length) {
        batch.size <- ceiling(nrow(x)/max.length)
        x <- aggregate(ts(x, frequency = batch.size), nfreq = 1, 
            FUN = mean)
    }
    else {
        batch.size <- 1
    }
    out <- do.spectrum0(x, max.freq = max.freq, order = order)
    out$spec <- out$spec * batch.size
    return(out)
}


do.spectrum0 <- function (x, max.freq = 0.5, order = 1) 
{
    fmla <- switch(order + 1, spec ~ one, spec ~ f1, spec ~ f1 + 
        f2)
    if (is.null(fmla)) 
        stop("invalid order")
    N <- nrow(x)
    Nfreq <- floor(N/2)
    freq <- seq(from = 1/N, by = 1/N, length = Nfreq)
    f1 <- sqrt(3) * (4 * freq - 1)
    f2 <- sqrt(5) * (24 * freq^2 - 12 * freq + 1)
    v0 <- numeric(ncol(x))
    for (i in 1:ncol(x)) {
        y <- x[, i]
        if (var(y) == 0) {
            v0[i] <- 0
        }
        else {
            yfft <- fft(y)
            spec <- Re(yfft * Conj(yfft))/N
            spec.data <- data.frame(one = rep(1, Nfreq), f1 = f1, 
                f2 = f2, spec = spec[1 + (1:Nfreq)], inset = I(freq <= 
                  max.freq))
            glm.out <- glm(fmla, family = Gamma(link = "log"), 
                data = spec.data, subset = spec.data$inset)
            v0[i] <- predict(glm.out, type = "response", newdata = data.frame(spec = 0, 
                one = 1, f1 = -sqrt(3), f2 = sqrt(5)))
        }
    }
    return(list(spec = v0))
}
