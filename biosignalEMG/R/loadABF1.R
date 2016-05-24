loadABF1 <- function(filename, const = 0.00489615) {
    if (missing(filename)) 
        stop("filename is not specified")
    if (!file.exists(filename)) 
        stop("Not such file or directory")
    if (file.access(filename, mode = 4) < 0) 
        stop("file is not readable")
    
    f <- file(filename, open = "rb")
    seek(f, 122)
    r <- readBin(f, what = "numeric", n = 1L, size = 4)
    seek(f, 10)
    L <- readBin(f, what = "integer", n = 1L, size = 4)
    seek(f, 120)
    nch <- readBin(f, what = "integer", n = 1L, size = 1)
    frec <- 1e+06/(nch * r)
    L <- (L%/%nch) * nch
    seek(f, 410)
    sseq <- readBin(f, what = "integer", n = 16, size = 2)
    seek(f, 442)
    snames <- readChar(f, nchars = rep(10, 16))
    seek(f, 602)
    sunits <- readChar(f, nchars = rep(8, 16))
    seek(f, 2048)
    data <- readBin(f, what = "integer", n = L, size = 2)
    data <- matrix(data * const, ncol = nch)
    close(f)
    
    return(emg(data, samplingrate = frec, units = sunits[sseq[1:nch] + 1], data.name = snames[sseq[1:nch] + 
        1]))
} 
