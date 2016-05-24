#'@title Bandwidth selection
#'
#'@description This function at fixk 0 finds the bandwidth, or lag length, of the
#' short, long, and run variance of a residual.
#'
#'@usage nw(v,fixk)
#'
#'@param v A Vector of error terms from a factor model
#'
#'@param fixk If fixk is 0, then this function will perform automatic bandwidth
#' selection. Otherwise, the integer placed here will be the selected bandwidth.
#'
#'@return k An integer that is the bandwidth chosen by \code{nw()}
#'
#'@return w The vector of long run variance
#'
#'@references Moon, R. & B. Perron (2004) Testing for a unit root in panels with
#' dynamic factors. Journal of Econometrics 122, 81-126.
#'

nw <- function(v, fixk) {
    
    Tn <- dim(v)[1]
    
    nreg <- dim(v)[2]
    
    rho <- NULL
    
    sigma <- NULL
    
    if (fixk == 0) {
        # auto bandwidth selection
        
        bot <- 0
        
        top <- 0
        
        for (i in 1:nreg) {
            rho[i] <- qr.solve(v[1:I(Tn - 1), i], v[2:Tn, i])
            
            e <- v[2:Tn, i] - rho[i] * v[1:I(Tn - 1), i]
            
            sigma[i] <- crossprod(e)/(Tn - 1)
            
            top <- top + 4 * (rho[i]^2) * (sigma[i]^2)/(((1 - rho[i])^6) * (1 + rho[i])^2)
            
            bot <- bot + (sigma[i]^2)/((1 - rho[i])^4)
        }
        
        alpha <- top/bot
        
        k <- ceiling(1.1447 * (alpha * Tn)^(1/3))
        # Trying Something
        if (k > I(Tn)) {
            k = Tn-2
        }
    } else {
        
        k <- fixk
    }
    
    w <- matrix(0, k, 1)
    
    for (i in 1:k) {
        
        x <- i/k
        
        w[i] <- 1 - i/(k + 1)
    }
    
    
    output <- list(k = k, w = w)
    return(output)
} 
