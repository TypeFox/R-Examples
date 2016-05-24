circFBM<- function(n, H, plotfBm=FALSE){
        if(missing(n)) n <- 500
        if(missing(H))
                H <- 0.6
        if(missing(plotfBm)) plotfBm <- 1       
## -------------------------------------------------------------------
## first line of the circulant matrix, C, built via covariances of fGn 
## -------------------------------------------------------------------
        lineC <- function(n, H, m){
        	k <- 0:(m - 1)
                H2 <- 2 * H
                v <- (abs((k - 1)/n)^H2 - 2 * (k/n)^H2 + ((k + 1)/n)^H2)/2
                ind <- c(0:(m/2 - 1), m/2, (m/2 - 1):1)
                v <- v[ind + 1]
                drop(v)
        }
##
## ---------------------
## next power of two > n
## ---------------------
        m <- 2
        repeat {
                m <- 2 * m
                if(m >= (n - 1))
                        break
        }
        stockm <- m     ##
## ----------------------------------------------
## research of the power of two (<2^18) such that
## C is definite positive
## ----------------------------------------------
        repeat {
                m <- 2 * m
                eigenvalC <- lineC(n, H, m)
                if (m<100) print(eigenvalC)
                eigenvalC <- Re(fft(c(eigenvalC), inverse = F))
                if((all(eigenvalC > 0)) | (m > 2^18))
                        break
        }
        if(m > 2^18) {
                cat("----> exact method, impossible!!", fill = T)
                cat("----> can't find m such that C is definite positive", fill
                         = T)
                break
        }
        else {
##
## -----------------------------------------------
## simulation of W=(Q)^t Z, where Z leads N(0,I_m)
## and  (Q)_{jk} = m^(-1/2) exp(-2i pi jk/m)
## -----------------------------------------------
                ar <- rnorm(m/2 + 1)
                ai <- rnorm(m/2 + 1)
                ar[1] <- sqrt(2) * ar[1]
                ar[(m/2 + 1)] <- sqrt(2) * ar[(m/2 + 1)]
                ai[1] <- 0
                ai[(m/2 + 1)] <- 0
                ar <- c(ar[c(1:(m/2 + 1))], ar[c((m/2):2)])
                aic <-  - ai
                ai <- c(ai[c(1:(m/2 + 1))], aic[c((m/2):2)])
                W <- complex(real = ar, imaginary = ai) ##
## -------------------------
## reconstruction of the fGn
## -------------------------
                W <- (sqrt(eigenvalC)) * W
                fGn <- fft(W, inverse = F)
                fGn <- (1/(sqrt(2 * m))) * fGn
                fGn <- Re(fGn[c(1:n)])
                fBm <- cumsum(fGn)
                fBm[1] <- 0     ##
## -----------
## plot of fBm
## -----------
                if(plotfBm) {
                        par(mfrow = c(1, 1))
                        time <- (0:(n - 1))/n
                        Nchar <- as.character(n)
                        Nleg <- paste(c("N= ", Nchar), collapse = " ")
                        Hchar <- as.character(round(H, 3))
                        Hleg <- paste(c(", H=", Hchar), collapse = "")
                        NHleg <- paste(c(Nleg, Hleg), collapse = "")
                        leg <- paste(c(
                                "Path of a fractional Brownian motion ----- parameters",
                                NHleg), collapse = " : ")
                        plot(time, fBm, type = "l", main = leg)
                }
                fBm
        }
}