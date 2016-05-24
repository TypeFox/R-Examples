###############################################################################
## finite-sampel under-/overshoot risk
###############################################################################

# n-fold convolution of Q(dx, dy)
conv.Q <- function(z, b, tau, K, n, m){
    if(z <= -n*b) return(0)
    if(z >= n*b) return(1)

    integrand.st <- function(x, b, tau){
        return(pnorm(b/abs(x) + tau*abs(x)) - pnorm(-b/abs(x) + tau*abs(x)))
    }
    stand <- E(K, integrand.st, b = b, tau = tau)
    
    if(n == 1){
        integrand.n1 <- function(x, z, b, tau){
            return(pnorm(z/abs(x) + tau*abs(x)) - pnorm(-b/abs(x) + tau*abs(x)))
        }
        erg <- E(K, integrand.n1, z = z, b = b, tau = tau)/stand

        return(erg)
    }
    
    M <- 2^m
    h <- 2*b/M
    z1 <- seq(from = -b, to = b, by = h)
    prob.x <- function(z1, tau, b, K){
        integrand.px <- function(x, z1, b, tau){
            return(pnorm(z1/abs(x) + tau*abs(x)) - pnorm(-b/abs(x) + tau*abs(x)))
        }
        erg <- E(K, integrand.px, z1 = z1, b = b, tau = tau)

        return(erg)                
    }
    p1 <- sapply(z1, prob.x, tau=tau, b=b, K=K)/stand
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
    x <- c(-n*b, seq(from = -n*b+n/2*h, to = n*b-n/2*h, by=h), n*b)
    pnfun1 <- approxfun(x = x+0.5*h, y = pn, yleft = 0, yright = pn[i.max+1])
    pnfun2 <- function(x) pnfun1(x) / pn[i.max+1]

    return(pnfun2(z))
}

setMethod("getFiRiskRegTS", signature(risk = "fiUnOvShoot",
                                      ErrorDistr = "Norm",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "ContNeighborhood"),
    function(risk, ErrorDistr, Regressor, neighbor, clip, stand, sampleSize,
             Algo, cont){
        eps <- neighbor@radius
        tau <- risk@width
        n <- sampleSize
        m <- distr::DefaultNrFFTGridPointsExponent

        if(Algo != "A"){
            if(cont == "left"){
                integrand.d1 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) + tau*abs(x)) + pnorm(-b/abs(x) - tau*abs(x)))
                }
                delta1 <- E(Regressor, integrand.d1, b = clip, tau = tau)
                delta1 <- (1-eps)*delta1 + eps
                K1 <- dbinom(0:n, size = n, prob = delta1)

                integrand.P1 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) - tau*abs(x)))
                }
                P1 <- E(Regressor, integrand.P1, b = clip, tau = tau)
                P1 <- (1-eps)*P1 + eps
                p1 <- P1/delta1

                summe1 <- numeric(n+1)
                summe1[1] <- 1 - conv.Q(z = 0, b = clip, tau = tau, K = Regressor, n = n, m = m)
                summe1[n+1] <- 1 - 0.5*(pbinom(q = n/2, size = n, prob = p1) 
                                        + pbinom(q = n/2-0.1, size = n, prob = p1))
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P1.ste <- sapply(z, conv.Q, b = clip, tau = tau, K = Regressor, n = n-k, m = m)
                    summe1[k+1] <- sum((1-P1.ste)*dbinom(j, size=k, prob=p1))
                }
                erg <- sum(summe1*K1)
            }else{
                integrand.d2 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) + tau*abs(x)) + pnorm(-b/abs(x) - tau*abs(x)))
                }
                delta2 <- E(Regressor, integrand.d2, b = clip, tau = tau)
                delta2 <- (1-eps)*delta2 + eps
                K2 <- dbinom(0:n, size = n, prob = delta2)

                integrand.P2 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) + tau*abs(x)))
                }
                P2 <- E(Regressor, integrand.P2, b = clip, tau = tau)
                P2 <- (1-eps)*P2
                p2 <- P2/delta2

                summe2 <- numeric(n+1)
                summe2[1] <- conv.Q(z = 0, b = clip, tau = -tau, K = Regressor, n = n, m = m)
                summe2[n+1] <- 0.5*(pbinom(q = n/2, size = n, prob = p2) 
                                    + pbinom(q = n/2-0.1, size = n, prob = p2))
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P2.ste <- sapply(z, conv.Q, b = clip, tau = -tau, K = Regressor, n = n-k, m = m)
                    summe2[k+1] <- sum(P2.ste*dbinom(j, size=k, prob=p2))
                }
                erg <- sum(summe2*K2)
            }
        }else{
            M <- 2^m
            h <- 2*clip/M
            z <- seq(from = -clip, to = clip, by = h)
        
            if(cont == "right"){
                integrand.px <- function(x, z, tau){pnorm(z/abs(x) + tau*abs(x))}
                prob.x <- function(z, tau, K, Ipx){
                    return(E(K, Ipx, z = z, tau = tau))
                }
                p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px)
                p1 <- (1-eps)*(p1[2:(M + 1)] - p1[1:M])
                p11 <- E(Regressor, integrand.px, z = -clip, tau = tau)
                p1[1] <- p1[1] + (1-eps)*p11
                p1M <- E(Regressor, integrand.px, z = -clip, tau = -tau)
                p1[M] <- p1[M] + (1-eps)*p1M + eps
            }else{
                integrand.px <- function(x, z, tau){pnorm(z/abs(x) - tau*abs(x))}
                prob.x <- function(z, tau, K, Ipx){
                    return(E(K, Ipx, z = z, tau = tau))
                }
                p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px)
                p1 <- (1-eps)*(p1[2:(M + 1)] - p1[1:M])
                p11 <- E(Regressor, integrand.px, z = -clip, tau = tau)
                p1[1] <- p1[1] + (1-eps)*p11 + eps
                p1M <- E(Regressor, integrand.px, z = -clip, tau = -tau)
                p1[M] <- p1[M] + (1-eps)*p1M
            }
        
            ## FFT
            pn <- c(p1, numeric((n-1)*M))
        
            ## convolution theorem for DFTs
            pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
            pn <- (abs(pn) >= .Machine$double.eps)*pn
            pn <- cumsum(pn)
        
            k <- n*(M-1)/2
            erg <- ifelse(n%%2 == 0, (pn[k]+pn[k+1])/2, pn[k+1])
            if(cont=="right") erg <- 1-erg
        }

        gc()
        return(list(fiUnOvShoot = erg))
    })

setMethod("getFiRiskRegTS", signature(risk = "fiUnOvShoot",
                                      ErrorDistr = "Norm",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "TotalVarNeighborhood"),
    function(risk, ErrorDistr, Regressor, neighbor, clip, stand, sampleSize,
             Algo, cont){
        delta <- neighbor@radius
        tau <- risk@width
        n <- sampleSize
        m <- distr::DefaultNrFFTGridPointsExponent

        if(Algo != "A"){
            if(cont == "left"){
                integrand.d11 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) - tau*abs(x)))
                }
                integrand.d12 <- function(x, b, tau){
                    return(pnorm(b/abs(x) - tau*abs(x)))
                }
                delta11 <- min(E(Regressor, integrand.d11, b = clip, tau = tau) + delta, 1)
                delta12 <- 1 - min(E(Regressor, integrand.d12, b = clip, tau = tau) + delta, 1)
                delta1 <- delta11 + delta12
                K1 <- dbinom(0:n, size=n, prob=delta1)

                integrand.P1 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) - tau*abs(x)))
                }
                P1 <- min(E(Regressor, integrand.P1, b = clip, tau = tau) + delta, 1)
                p1 <- min(P1/delta1, 1)

                summe1 <- numeric(n+1)
                summe1[1] <- 1 - conv.Q(z = 0, b = clip, tau = tau, K = Regressor, n = n, m = m)
                summe1[n+1] <- 1 - 0.5*(pbinom(q = n/2, size = n, prob = p1)
                                        + pbinom(q = n/2-0.1, size = n, prob = p1))
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P1.ste <- sapply(z, conv.Q, b = clip, tau = tau, K = Regressor, n = n-k, m = m)
                    summe1[k+1] <- sum((1-P1.ste)*dbinom(j, size = k, prob = p1))
                }
                erg <- sum(summe1*K1)
            }else{
                integrand.d21 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) + tau*abs(x)))
                }
                integrand.d22 <- function(x, b, tau){
                    return(pnorm(b/abs(x) + tau*abs(x)))
                }
                delta21 <- max(0, E(Regressor, integrand.d21, b = clip, tau = tau) - delta)
                delta22 <- 1 - max(0, E(Regressor, integrand.d22, b = clip, tau = tau) - delta)
                delta2 <- delta21 + delta22
                K2 <- dbinom(0:n, size = n, prob = delta2)

                integrand.P2 <- function(x, b, tau){
                    return(pnorm(-b/abs(x) + tau*abs(x)))
                }
                P2 <- max(0, E(Regressor, integrand.P2, b = clip, tau = tau) - delta)
                p2 <- P2/delta2

                summe2 <- numeric(n+1)
                summe2[1] <- conv.Q(z = 0, b = clip, tau = -tau, K = Regressor, n = n, m = m)
                summe2[n+1] <- 0.5*(pbinom(q = n/2, size = n, prob = p2) 
                                    + pbinom(q = n/2-0.1, size = n, prob = p2))
                for(k in 1:(n-1)){
                    j <- 0:k
                    z <- clip*(k-2*j)
                    P2.ste <- sapply(z, conv.Q, b = clip, tau = -tau, K = Regressor, n = n-k, m = m)
                    summe2[k+1] <- sum(P2.ste*dbinom(j, size = k, prob = p2))
                }
                erg <- sum(summe2*K2)
            }
        }else{
            M <- 2^m
            h <- 2*clip/M
            z <- seq(from = -clip, to = clip, by = h)
        
            if(cont == "right"){
                integrand.px <- function(x, z, tau){pnorm(z/abs(x) + tau*abs(x))}
                prob.x <- function(z, tau, K, Ipx){
                    return(E(K, Ipx, z = z, tau = tau))
                }
                p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px)
                p1 <- p1[2:(M + 1)] - p1[1:M]
                p11 <- E(Regressor, integrand.px, z = -clip, tau = tau)
                p1[1] <- p1[1] + p11 - delta
                p1M <- E(Regressor, integrand.px, z = -clip, tau = -tau)
                p1[M] <- p1[M] + p1M + delta
            }else{
                integrand.px <- function(x, z, tau){pnorm(z/abs(x) - tau*abs(x))}
                prob.x <- function(z, tau, K, Ipx){
                    return(E(K, Ipx, z = z, tau = tau))
                }
                p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px)
                p1 <- p1[2:(M + 1)] - p1[1:M]
                p11 <- E(Regressor, integrand.px, z = -clip, tau = tau)
                p1[1] <- p1[1] + p11 + delta
                p1M <- E(Regressor, integrand.px, z = -clip, tau = -tau)
                p1[M] <- p1[M] + p1M - delta
            }
        
            ## FFT
            pn <- c(p1, numeric((n-1)*M))
        
            ## convolution theorem for DFTs
            pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
            pn <- (abs(pn) >= .Machine$double.eps)*pn
            pn <- cumsum(pn)
        
            k <- n*(M-1)/2
            erg <- ifelse(n%%2 == 0, (pn[k]+pn[k+1])/2, pn[k+1])
            if(cont=="right") erg <- 1-erg
        }

        gc()
        return(list(fiUnOvShoot = erg))
    })

setMethod("getFiRiskRegTS", signature(risk = "fiUnOvShoot",
                                      ErrorDistr = "Norm",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "CondContNeighborhood"),
    function(risk, ErrorDistr, Regressor, neighbor, clip, stand, sampleSize, cont){
        eps <- neighbor@radiusCurve
        tau <- risk@width
        n <- sampleSize
        m <- distr::DefaultNrFFTGridPointsExponent

        if(is(Regressor, "AbscontDistribution")){
            xlower <- ifelse(is.finite(q(Regressor)(0)), q(Regressor)(0), q(Regressor)(distr::TruncQuantile))
            xupper <- ifelse(is.finite(q(Regressor)(1)), q(Regressor)(1), q(Regressor)(1 - distr::TruncQuantile))
            x.vec <- seq(from = xlower, to = xupper, length = 1000)
        }else{
            if(is(Regressor, "DiscreteDistribution"))
                x.vec <- support(Regressor) 
            else
                x.vec <- unique(r(Regressor)(distr::RtoDPQ.e))
        }

        b.vec <- sapply(x.vec, clip)
        Max <- max(b.vec*abs(x.vec), na.rm = TRUE)
        M <- 2^m
        h <- 2*Max/M
        z <- seq(from = -Max, to = Max, by = h)
        
        if(cont == "right"){
            integrand.px <- function(x, z, tau, b, eps){
                if(z < -abs(x)*b(x)) return(0)
                if(z == -abs(x)*b(x)) return((1-eps(x))*pnorm(-b(x) + tau*abs(x)))
                if(z == abs(x)*b(x)) return((1-eps(x))*pnorm(-b(x) - tau*abs(x)) + eps(x))
                if(z > abs(x)*b(x)) return(1)
                return((1-eps(x))*pnorm(z/abs(x) + tau*abs(x)))
            }
            prob.x <- function(z, tau, K, Ipx, b, eps){
                return(E(K, Ipx, z = z, tau = tau, b = b, eps = eps))
            }
            p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px, b = clip, eps = eps)
            p1 <- c(p1[2], p1[3:M] - p1[2:(M-1)])
            p1 <- c(p1, 1-sum(p1))
        }else{
            integrand.px <- function(x, z, tau, b, eps){
                if(z < -abs(x)*b(x)) return(0)
                if(z == -abs(x)*b(x)) return((1-eps(x))*pnorm(-b(x) - tau*abs(x)) + eps(x))
                if(z == abs(x)*b(x)) return((1-eps(x))*pnorm(-b(x) + tau*abs(x)))
                if(z > abs(x)*b(x)) return(1)
                return((1-eps(x))*pnorm(z/abs(x) - tau*abs(x)) + eps(x))
            }
            prob.x <- function(z, tau, K, Ipx, b, eps){
                return(E(K, Ipx, z = z, tau = tau, b = b, eps = eps))
            }
            p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px, b = clip, eps = eps)
            p1 <- c(p1[2], p1[3:M] - p1[2:(M-1)])
            p1 <- c(p1, 1-sum(p1))
        }
        
        ## FFT
        pn <- c(p1, numeric((n-1)*M))
        
        ## convolution theorem for DFTs
        pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
        pn <- (abs(pn) >= .Machine$double.eps)*pn
        pn <- cumsum(pn)
        
        k <- n*(M-1)/2
        erg <- ifelse(n%%2 == 0, (pn[k]+pn[k+1])/2, pn[k+1])
        if(cont=="right") erg <- 1-erg
        

        gc()
        return(list(fiUnOvShoot = erg))
    })

setMethod("getFiRiskRegTS", signature(risk = "fiUnOvShoot",
                                      ErrorDistr = "Norm",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "CondTotalVarNeighborhood"),
    function(risk, ErrorDistr, Regressor, neighbor, clip, stand, sampleSize, cont){
        delta <- neighbor@radiusCurve
        tau <- risk@width
        n <- sampleSize
        m <- distr::DefaultNrFFTGridPointsExponent

        if(is(Regressor, "AbscontDistribution")){
            xlower <- ifelse(is.finite(q(Regressor)(0)), q(Regressor)(0), q(Regressor)(distr::TruncQuantile))
            xupper <- ifelse(is.finite(q(Regressor)(1)), q(Regressor)(1), q(Regressor)(1 - distr::TruncQuantile))
            x.vec <- seq(from = xlower, to = xupper, length = 1000)
        }else{
            if(is(Regressor, "DiscreteDistribution"))
                x.vec <- support(Regressor) 
            else
                x.vec <- unique(r(Regressor)(distr::RtoDPQ.e))
        }

        b.vec <- sapply(x.vec, clip)
        Max <- max(b.vec*abs(x.vec), na.rm = TRUE)
        M <- 2^m
        h <- 2*Max/M
        z <- seq(from = -Max, to = Max, by = h)
        
        if(cont == "right"){
            integrand.px <- function(x, z, tau, b, delta){
                if(z < -abs(x)*b(x)) return(0)
                if(z == -abs(x)*b(x)) return(max(0, pnorm(-b(x) + tau*abs(x)) - delta(x)))
                if(z == abs(x)*b(x)) return(1 - max(0, pnorm(b(x) + tau*abs(x)) - delta(x)))
                if(z > abs(x)*b(x)) return(1)
                return(max(0, pnorm(z/abs(x) + tau*abs(x)) - delta(x)))
            }
            prob.x <- function(z, tau, K, Ipx, b, delta){
                return(E(K, Ipx, z = z, tau = tau, b = b, delta = delta))
            }
            p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px, b = clip, delta = delta)
            p1 <- c(p1[2], p1[3:M] - p1[2:(M-1)])
            p1 <- c(p1, 1-sum(p1))
        }else{
            integrand.px <- function(x, z, tau, b, delta){
                if(z < -abs(x)*b(x)) return(0)
                if(z == -abs(x)*b(x)) return(min(1, pnorm(-b(x) - tau*abs(x)) + delta(x)))
                if(z == abs(x)*b(x)) return(1 - min(1, pnorm(b(x) - tau*abs(x)) + delta(x)))
                if(z > abs(x)*b(x)) return(1)
                return(min(1, pnorm(z/abs(x) - tau*abs(x)) + delta(x)))
            }
            prob.x <- function(z, tau, K, Ipx, b, delta){
                return(E(K, Ipx, z = z, tau = tau, b = b, delta = delta))
            }
            p1 <- sapply(z, prob.x, tau = tau, K = Regressor, Ipx = integrand.px, b = clip, delta = delta)
            p1 <- c(p1[2], p1[3:M] - p1[2:(M-1)])
            p1 <- c(p1, 1-sum(p1))
        }
        
        ## FFT
        pn <- c(p1, numeric((n-1)*M))
        
        ## convolution theorem for DFTs
        pn <- Re(fft(fft(pn)^n, inverse = TRUE)) / (n*M)
        pn <- (abs(pn) >= .Machine$double.eps)*pn
        pn <- cumsum(pn)
        
        k <- n*(M-1)/2
        erg <- ifelse(n%%2 == 0, (pn[k]+pn[k+1])/2, pn[k+1])
        if(cont=="right") erg <- 1-erg
        
        gc()
        return(list(fiUnOvShoot = erg))
    })
