syntheticemg <- function(n.length.out = 10000, on.sd = 1, on.duration.mean = 350, 
    on.duration.sd = 10, off.sd = 0.05, off.duration.mean = 300, off.duration.sd = 20, 
    on.mode.pos = 0.75, shape.factor = 0.5, samplingrate = 0, units = "", data.name = "Synthetic EMG") {
    n.length.out <- round(n.length.out)
    if (n.length.out <= 0) 
        stop("length must be non-negative number")
    if (on.sd <= 0) 
        stop("on.sd (standard deviation of the active phases) must be non-negative number")
    if (off.sd <= 0) 
        stop("off.sd (standard deviation of the inactive phases) must be non-negative number")
    if (off.sd >= on.sd) 
        stop("on.sd (standard deviation of the active phases) must be greater than off.sd (standard deviation of the inactive phases)")
    if (on.duration.sd <= 0) 
        stop("on.sd (standard deviation of the active phases) must be non-negative number")
    if (off.duration.sd <= 0) 
        stop("off.sd (standard deviation of the inactive phases) must be non-negative number")
    if (on.mode.pos < 0 | on.mode.pos > 1) 
        stop("'on.mode.pos' outside [0,1]")
    
    i <- 0
    b <- numeric()
    semg <- numeric()
    while (i < n.length.out) {
        n.off <- round(rnorm(1, off.duration.mean, off.duration.sd))
        if (n.off < 1) 
            n.off <- 1
        n.on <- round(rnorm(1, on.duration.mean, on.duration.sd))
        if (n.on < 1) 
            n.on <- 1
        b <- c(b, rep(0, n.off), rep(1, n.on))
        semg <- c(semg, rnorm(n.off, 0, off.sd), c(seq(off.sd/on.sd, 1, length.out = round(on.mode.pos * 
            n.on)), seq(1, off.sd/on.sd, length.out = n.on - round(on.mode.pos * 
            n.on)))^(shape.factor) * rnorm(n.on, 0, on.sd))
        i <- i + n.off + n.on
    }
    b <- head(b, n.length.out)
    semg <- head(semg, n.length.out)
    object <- list(values = semg, units = units, samplingrate = samplingrate, data.name = data.name, 
        on.off = b)
    class(object) <- "emg"
    return(object)
} 
