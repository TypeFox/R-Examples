###############################################################################
## asymptotic under-/overshoot risk
###############################################################################
setMethod("getRiskIC", signature(IC = "TotalVarIC",
                                 risk = "asUnOvShoot",
                                 neighbor = "UncondNeighborhood",
                                 L2Fam = "missing"),
    function(IC, risk, neighbor){
        radius <- neighbor@radius
        L2Fam <- eval(IC@CallL2Fam)
        L2deriv <- L2Fam@L2derivDistr[[1]]
        if((length(L2Fam@L2derivDistr) > 1) | !is(L2deriv, "UnivariateDistribution"))
            stop("restricted to 1-dimensional parameteric models")

        bound <- risk@width*(-m1df(L2deriv, 0))
        if(is(neighbor, "ContNeighborhood")){
            if(radius > 2*bound)
                stop("boundedness condition is violated!")
            if(radius == 2*bound){
                zi <- sign(as.vector(trafo))
                A <- as.matrix(zi)
                b <- zi*as.vector(trafo)*2*risk@width/radius
                p0 <- p(L2deriv)(0)
                if(is(L2deriv, "AbscontDistribution"))
                    ws0 <- 0
                else
                    ws0 <- d(L2deriv)(0)

                if(zi == 1)
                    a <- -b*(1-p0)/(1-ws0)
                else
                    a <- b*(p0-ws0)/(1-ws0)

                asCov <- a^2*(p0-ws0) + (zi*a+b)^2*(1-p0)
                erg <- pnorm(-risk@width*sqrt(asCov))
            }
        }

        if(is(neighbor, "TotalVarNeighborhood")){
            if(radius > bound)
                stop("boundedness condition is violated!")
            if(radius == bound){
                zi <- sign(as.vector(trafo))
                A <- as.matrix(zi)
                b <- zi*as.vector(trafo)*risk@width/radius
                p0 <- p(L2deriv)(0)
                if(is(L2deriv, "AbscontDistribution"))
                    ws0 <- 0
                else
                    ws0 <- d(L2deriv)(0)

                if(zi == 1)
                    a <- -b*(1-p0)/(1-ws0)
                else
                    a <- b*(p0-ws0)/(1-ws0)

                asCov <- a^2*(p0-ws0) + (zi*a+b)^2*(1-p0)
                erg <- pnorm(-risk@width*sqrt(asCov))
            }
        }

        stand <- as.vector(stand(IC))
        g0 <- clipLo(IC)/abs(stand)
        c0 <- clipUp(IC)/abs(stand) - g0
        s <- sqrt(g0^2*p(L2deriv)(g0)
                  + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
                  + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0))
        erg <- pnorm(-risk@width*s)

        nghb <- paste(neighbor@type, "with radius", neighbor@radius)

        return(list(asUnOvShoot = list(distribution = .getDistr(L2Fam),
                    neighborhood = nghb, value = erg)))
    })
###############################################################################
## finite-sample under-/overshoot risk
###############################################################################
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "fiUnOvShoot",
                                 neighbor = "ContNeighborhood",
                                 L2Fam = "missing"),
    function(IC, risk, neighbor, sampleSize, Algo = "A", cont = "left"){
        L2Fam <- eval(IC@CallL2Fam)
        Distr <- L2Fam@distribution
        if(!is(Distr, "Norm"))
            stop("restricted to 1-dimensional normal location")

        eps <- neighbor@radius
        tau <- risk@width

        if(!(is(IC, "ContIC") | is(IC, "TotalVarIC")))
            stop("'IC' has to be of class 'ContIC' or 'TotalVarIC'")
        if(is(IC, "ContIC"))
            clip <- clip(IC)/as.vector(stand(IC))
        if(is(IC, "TotalVarIC"))
            clip <- clipUp(IC)/as.vector(stand(IC))

        n <- sampleSize
        m <- getdistrOption("DefaultNrFFTGridPointsExponent")

        if(eps >= 1 - 1/(2*pnorm(risk@width))){
            warning("disjointness condition is violated!")
            erg <- 0.5
        }else{
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
        }

        nghb <- paste(neighbor@type, "with radius", neighbor@radius)

        return(list(fiUnOvShoot = list(distribution = .getDistr(eval(IC@CallL2Fam)),
                    neighborhood = nghb, value = erg)))
    })
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "fiUnOvShoot",
                                 neighbor = "TotalVarNeighborhood",
                                 L2Fam = "missing"),
    function(IC, risk, neighbor, sampleSize, Algo = "A", cont = "left"){
        L2Fam <- eval(IC@CallL2Fam)
        Distr <- L2Fam@distribution
        if(!is(Distr, "Norm"))
            stop("restricted to 1-dimensional normal location")

        delta <- neighbor@radius
        tau <- risk@width

        if(!(is(IC, "ContIC") | is(IC, "TotalVarIC")))
            stop("'IC' has to be of class 'ContIC' or 'TotalVarIC'")
        if(is(IC, "ContIC"))
            clip <- clip(IC)/as.vector(stand(IC))
        if(is(IC, "TotalVarIC"))
            clip <- clipUp(IC)/as.vector(stand(IC))

        n <- sampleSize
        m <- getdistrOption("DefaultNrFFTGridPointsExponent")

        if(delta >= pnorm(risk@width) - 0.5){
            warning("disjointness condition is violated!")
            erg <- 0.5
        }else{
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
        }

        nghb <- paste(neighbor@type, "with radius", neighbor@radius)

        return(list(fiUnOvShoot = list(distribution = .getDistr(eval(IC@CallL2Fam)),
                    neighborhood = nghb, value = erg)))
    })
