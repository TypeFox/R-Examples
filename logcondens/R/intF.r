intF <- function (s, res){

    x <- res$x
    phi <- res$phi
    Fhat <- res$Fhat    
    
    if (min(s) < min(x) | max(s) > max(x)){cat("All elements of s must be in [x_1, x_n]!")} else {
        n <- length(x)
        dx <- c(NA, diff(x))
        dphi <- c(NA, diff(phi))
        slope <- dphi / dx
        f <- exp(phi)
        intF.xi <- c(0, rep(NA, n - 1))

        for (i in 2:n){
            if (abs(slope[i]) <= 1e-6){inte <- f[i - 1] * dx[i] / 2}
            if (abs(slope[i]) > 1e-6){inte <- slope[i] ^ (-1) * (J00(phi[i - 1], phi[i], 1) - f[i - 1])}
            intF.xi[i] <- dx[i] * (Fhat[i - 1] + inte)
            }

        intF.xi <- cumsum(intF.xi)
        intF.s <- rep(NA, length(s))
        for (k in 1:length(s)){
            x0 <- s[k]
            j <- max((1:n)[x <= x0])
            j <- min(j, n - 1)
            xj <- x[j]
            
            if (abs(slope[j + 1]) <= 1e-6){
                  tmp <- f[j] * (x0 - xj) ^ 2 / 2}
            
            if (abs(slope[j + 1]) > 1e-6){
                  tmp <- dx[j + 1] * (slope[j + 1] ^ (-1) * J00(phi[j], 
                    phi[j + 1], (x0 - xj) / dx[j + 1]) - (x0 - xj) / dphi[j + 1] * f[j])
                    }
                  
            intF.s[k] <- intF.xi[j] + (x0 - xj) * Fhat[j] + tmp
        }
        
        return(intF.s)
    }
}

