###############################################################################
## finite-sample under-/overshoot risk
###############################################################################

# cdf of truncated normal distribution
ptnorm <- function(x, mu, A, B){
    ((A <= x)*(x <= B)*(pnorm(x-mu)-pnorm(A-mu))/(pnorm(B-mu)-pnorm(A-mu))
    + (x > B))
}

# n-fold convolution for truncated normal distributions
conv.tnorm <- function(z, A, B, mu, n, m){
    if(n == 1) return(ptnorm(z, mu = mu, A = A, B = B))
    if(z <= n*A) return(0)
    if(z >= n*B) return(1)
    
    M <- 2^m
    h <- (B-A)/M
    x <- seq(from = A, to = B, by = h)
    p1 <- ptnorm(x, mu = mu, A = A, B = B)
    p1 <- p1[2:(M + 1)] - p1[1:M]

    ## FFT
    pn <- c(p1, numeric((n-1)*M))

    ## convolution theorem for DFTs
    pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
    pn <- (abs(pn) >= .Machine$double.eps)*pn
    i.max <- n*M-(n-2)
    pn <- c(0,pn[1:i.max])
    pn <- cumsum(pn)

    ## cdf with continuity correction h/2
    x <- c(n*A,seq(from = n*A+n/2*h, to = n*B-n/2*h, by=h),n*B)
    pnfun1 <- approxfun(x = x+0.5*h, y = pn, yleft = 0, yright = pn[i.max+1])
    pnfun2 <- function(x) pnfun1(x) / pn[i.max+1]

    return(pnfun2(z))
}


setMethod("getFiRisk", signature(risk = "fiUnOvShoot",
                                 Distr = "Norm",
                                 neighbor = "ContNeighborhood"),
    function(risk, Distr, neighbor, clip, stand, sampleSize, Algo, cont){
        eps <- neighbor@radius
        tau <- risk@width
        n <- sampleSize
        m <- getdistrOption("DefaultNrFFTGridPointsExponent")
        
        if(Algo == "B"){
            if(cont == "left"){
                delta1 <- (1-eps)*(pnorm(-clip+tau) + pnorm(-clip-tau)) + eps
                K1 <- dbinom(0:n, size = n, prob = delta1)
                P1 <- (1-eps)*pnorm(-clip-tau) + eps
                p1 <- P1/delta1

                summe1 <- numeric(n+1)
                summe1[1] <- 1 - conv.tnorm(z = 0, A = -clip, B = clip, mu = -tau, n = n, m = m)
                summe1[n+1] <- (1 - 0.5*(pbinom(q = n/2, size = n, prob = p1) 
                                + pbinom(q = n/2-0.1, size = n, prob = p1)))
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P1.ste <- sapply(z, conv.tnorm, A = -clip, B = clip, mu = -tau, n = n-k, m = m)
                    summe1[k+1] <- sum((1-P1.ste)*dbinom(j, size = k, prob = p1))
                }
                erg <- sum(summe1*K1)
            }else{
                delta2 <- (1-eps)*(pnorm(-clip+tau) + pnorm(-clip-tau)) + eps
                K2 <- dbinom(0:n, size = n, prob = delta2)
                P2 <- (1-eps)*pnorm(-clip+tau)
                p2 <- P2/delta2

                summe2 <- numeric(n+1)
                summe2[1] <- conv.tnorm(z = 0, A = -clip, B = clip, mu = tau, n = n, m = m)
                summe2[n+1] <- 0.5*(pbinom(q = n/2, size = n, prob = p2) 
                                    + pbinom(q = n/2-0.1, size = n, prob = p2))
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P2.ste <- sapply(z, conv.tnorm, A = -clip, B = clip, mu = tau, n = n-k, m = m)
                    summe2[k+1] <- sum(P2.ste*dbinom(j, size=k, prob=p2))
               }
                erg <- sum(summe2*K2)
            }
        }else{
            M <- 2^m
            h <- 2*clip/M
            x <- seq(from = -clip, to = clip, by = h)

            if(cont == "right"){
                p1 <- pnorm(x+tau)
                p1 <- (1-eps)*(p1[2:(M + 1)] - p1[1:M])
                p1[1] <- p1[1] + (1-eps)*pnorm(-clip+tau)
                p1[M] <- p1[M] + (1-eps)*pnorm(-clip-tau) + eps
            }else{
                p1 <- pnorm(x-tau)
                p1 <- (1-eps)*(p1[2:(M + 1)] - p1[1:M])
                p1[1] <- p1[1] + (1-eps)*pnorm(-clip-tau) + eps
                p1[M] <- p1[M] + (1-eps)*pnorm(-clip+tau)
            }
        
            ## FFT
            pn <- c(p1, numeric((n-1)*M))

            ## convolution theorem for DFTs
            pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
            pn <- (abs(pn) >= .Machine$double.eps)*pn
            pn <- cumsum(pn)

            k <- n*(M-1)/2
            erg <- ifelse(n%%2 == 0, (pn[k]+pn[k+1])/2, pn[k+1])
            if(cont == "right") erg <- 1 - erg
        }

        return(list(fiUnOvShoot = erg))
    })

setMethod("getFiRisk", signature(risk = "fiUnOvShoot",
                                 Distr = "Norm",
                                 neighbor = "TotalVarNeighborhood"),
    function(risk, Distr, neighbor, clip, stand, sampleSize, Algo, cont){
        delta <- neighbor@radius
        tau <- risk@width
        n <- sampleSize
        m <- getdistrOption("DefaultNrFFTGridPointsExponent")

        if(Algo == "B"){
            if(cont == "left"){
                delta1 <- min(pnorm(-clip-tau)+delta, 1) + 1 - min(pnorm(clip-tau)+delta, 1)
                K1 <- dbinom(0:n, size = n, prob = delta1)
                P1 <- min(pnorm(-clip-tau) + delta, 1)
                p1 <- min(P1/delta1, 1)

                summe1 <- numeric(n+1)
                summe1[1] <- 1 - conv.tnorm(z = 0, A = -clip, B = clip, mu = -tau, n = n, m = m)
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P1.ste <- sapply(z, conv.tnorm, A = -clip, B = clip, mu = -tau, n = n-k, m = m)
                    summe1[k+1] <- sum((1-P1.ste)*dbinom(j, size = k, prob = p1))
                }
                summe1[n+1] <- 1 - 0.5*(pbinom(q = n/2, size = n, prob = p1)
                                        + pbinom(q = n/2-0.1, size = n, prob = p1))
                erg <- sum(summe1*K1)
            }else{
                delta2 <- max(0, pnorm(-clip+tau)-delta) + 1 - max(0, pnorm(clip+tau)-delta)
                K2 <- dbinom(0:n, size = n, prob = delta2)
                P2 <- max(0, pnorm(-clip+tau) - delta)
                p2 <- P2/delta2

                summe2 <- numeric(n+1)
                summe2[1] <- conv.tnorm(z = 0, A = -clip, B = clip, mu = tau, n = n, m = m)
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P2.ste <- sapply(z, conv.tnorm, A = -clip, B = clip, mu = tau, n = n-k, m = m)
                    summe2[k+1] <- sum(P2.ste*dbinom(j, size = k, prob = p2))
                }
                summe2[n+1] <- 0.5*(pbinom(q = n/2, size = n, prob = p2) 
                                    + pbinom(q = n/2-0.1, size = n, prob = p2))
                erg <- sum(summe2*K2)
            }
        }else{
            M <- 2^m
            h <- 2*clip/M
            x <- seq(from = -clip, to = clip, by = h)

            if(cont == "right"){
                p1 <- pnorm(x+tau)
                p1 <- p1[2:(M + 1)] - p1[1:M]
                p1[1] <- p1[1] + pnorm(-clip+tau) - delta
                p1[M] <- p1[M] + pnorm(-clip-tau) + delta
            }else{
                p1 <- pnorm(x-tau)
                p1 <- p1[2:(M + 1)] - p1[1:M]
                p1[1] <- p1[1] + pnorm(-clip-tau) + delta
                p1[M] <- p1[M] + pnorm(-clip+tau) - delta
            }

            ## FFT
            pn <- c(p1, numeric((n-1)*M))

            ## convolution theorem for DFTs
            pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
            pn <- (abs(pn) >= .Machine$double.eps)*pn
            pn <- cumsum(pn)
    
            k <- n*(M-1)/2
            erg <- ifelse(n%%2 == 0, (pn[k]+pn[k+1])/2, pn[k+1])
            if(cont == "right") erg <- 1-erg
        }

        return(list(fiUnOvShoot = erg))
    })
