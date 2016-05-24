onoff_bonato <- function(data, channel, sigma_n, Pfa = 0.05, m = 5, r0 = 1, minL = 15, 
    data.name) {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.emg(data)) 
        stop("an object of class 'emg' is required")
    if (missing(sigma_n)) 
        stop("An estimation of the standard deviation 'sigma_n' of the baseline noise is requiered")
    if (missing(channel)) {
        if (missing(data.name)) 
            data <- extractchannel(data) else data <- extractchannel(data, data.name = data.name)
    } else {
        if (missing(data.name)) 
            data <- extractchannel(data, channel) else data <- extractchannel(data, channel, data.name)
    }
    if (is.null(m)) 
        m <- 5  #(for 10ms)
    if (is.null(minL)) 
        minL <- 15  #(for 30ms)
    
    Pzeta <- 1e-04
    while (1 - pbinom(r0 - 1, m, Pzeta) < Pfa) {
        Pzeta <- Pzeta + 1e-04
    }
    zeta <- (-2 * sigma_n^2 * log(Pzeta))
    
    z <- tail(data$values, 1)^2 + head(data$values, -1)^2
    cs <- cumsum(z > zeta)
    pcs <- tail(cs, -(m - 1)) - c(0, head(cs, -m))
    detected <- c(pcs >= r0, rep(0, m))
    #return(detected)
    cad <- rle(detected)
    
    cadL <- cad$lengths[1]
    cadV <- cad$values[1]
    k <- 1
    i <- 2
    L <- length(cad$values)
    while (i < L) {
        if (cad$lengths[i] < minL) {
            cadL[k] <- cadL[k] + cad$lengths[i] + cad$lengths[i + 1]
            i <- i + 2
        } else {
            cadL[k + 1] <- cad$lengths[i]
            cadV[k + 1] <- cad$values[i]
            k <- k + 1
            i <- i + 1
        }
    }
    if (i < L) 
        i <- i + 1
    if (i == L) {
        if (cad$lengths[i] < minL) {
            cadL[k] <- cadL[k] + cad$lengths[i]
        } else {
            cadL[k + 1] <- cad$lengths[i]
            cadV[k + 1] <- cad$values[i]
        }
    }
    if (cadL[1] < minL) {
        cadL[2] <- cadL[2] + cadL[1]
        cadL <- tail(cadL, -1)
        cadV <- tail(cadV, -1)
    }
    cad2 <- list(lengths = cadL, values = cadV)
    class(cad2) <- "rle"
    detected <- inverse.rle(cad2)
    return(detected)
} 
