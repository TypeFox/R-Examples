#'@title Finding optimal lag for dfgls test
#'
#'@description This function finds the optimal lag kstar for the dfgls
#' test.
#'
#'@usage s2ar(yts, penalty, kmax, kmin)
#'
#'@param yts A NxT matrix containing the data to find kstar for.
#'
#'@param penalty a binary selection of 0 or 1. 0 uses the MAIC, a penalty on
#' k that accounts for the bias in the sum of the autoregressive coefficients.
#'  1 uses the more general form MIC.
#'
#'@param  kmax The maximum number of lags for the vector autoregressions. An
#' upper bound of \code{(12x(T/100)^.25)^8} is suggested
#' in Schwert (1989)
#'
#'@param  kmin The minimum number of lags for the vector autoregression. k equal to 0
#' is a reasonable point.
#'
#'@return kstar A vector of optimal lags for each column of yts
#'
#'@references Schwert, G. W. 1989. Tests for unit roots: A Monte Carlo
#' investigation. Journal of Business and Economic Statistics 2: 147-159.
#'
#'            Serana Ng and P. Perron. 2000. Lag length selection and the construction of unit root
#' tests with good size and power. Econometrica 69:1519-1554.
#'
#'

s2ar <- function(yts, penalty, kmax, kmin) {
    
    nt <- nrow(yts)
    
    minn <- 9999999999
    
    tauu <- matrix(0, I(kmax + 1), 1)
    
    s2e <- 999 * matrix(1, I(kmax + 1), 1)
    
    dyts <- mydiff(yts, 1)
    
    reg <- lagn(yts, 1)
    
    for (i in 1:kmax) {
        
        reg <- cbind(reg, lagn(dyts, i))
    }
    
    dyts0 <- dyts
    
    reg0 <- reg
    
    dyts0 <- trimr(dyts, I(kmax + 1), 0)
    
    reg0 <- trimr(reg, I(kmax + 1), 0)
    
    sumy <- sum(reg0[, 1] * reg0[, 1])
    
    nef <- nt - kmax - 1
    for (k in kmin:kmax) {
        
        b <- myols(reg0[, 1:I(k + 1)], dyts0)
        
        e <- dyts0 - reg0[, 1:I(k + 1)] %*% b
        
        s2e[I(k + 1), ] <- t(e) %*% e/nef
        
        tauu[I(k + 1), ] <- (b[1] * b[1]) * sumy/s2e[I(k + 1)]
    }
    
    kk <- seq(0, kmax)
    
    if (penalty == 0) {
        
        mic <- log(s2e) + 2 * (kk + tauu)/nef
        
    } else {
        
        mic <- log(s2e) + log(nef) * (kk)/nef
        
    }
    
    kstar <- minindc(mic) - 1
    
    return(kstar)
} 
