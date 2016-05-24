PR <- function(b, bE, t) {
    if (!is.vector(b) | !is.vector(bE)) 
        stop("Function is only defined for vectors")
    if (length(b) != length(bE)) 
        stop("Function is only defined for vectors of the same length")
    if (length(b) < 1) 
        stop("Function is only defined for non-empty vectors")
    if (t <= 0) 
        stop("Tolerance value 't' should be positive")
    
    n <- length(b)
    CPs_b <- c(0, 2^b[2:n] - 2^b[1:(n - 1)])
    ps_b <- which(CPs_b != 0)
    Nevents <- length(ps_b)
    
    CPs <- c(0, 2^bE[2:n] - 2^bE[1:(n - 1)])
    ps <- which(CPs != 0)
    Ndet <- length(ps)
    
    if (Nevents & Ndet != 0) {
        NTP <- 0
        NFP <- 0
        for (i in 1:length(ps)) {
            k <- ps[i]
            tipo = 2^bE[k] - 2^bE[k - 1]
            ps <- which(sign(CPs_b) == sign(tipo))
            d <- min(abs(ps - k))
            if (d < t) 
                NTP <- NTP + 1 else NFP <- NFP + 1
        }
        TPR <- NTP/Nevents
        FPR <- NFP/Ndet
    } else {
        TPR <- 1
        FPR <- 0
        if (Nevents != 0) {
            TPR <- 0
            FPR <- 0
        }
        if (Ndet != 0) {
            TPR <- 0
            FPR <- 1
        }
        
    }
    return(list(TPR = TPR, FPR = FPR))
} 
