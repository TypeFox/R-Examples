###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getFixRobRegTypeIC", signature(ErrorDistr = "Norm",
                                          Regressor = "UnivariateDistribution", 
                                          risk = "fiUnOvShoot", 
                                          neighbor = "UncondNeighborhood"),
    function(ErrorDistr, Regressor, risk, neighbor, sampleSize, upper, maxiter, 
             tol, warn, Algo, cont){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE))
            stop("'radius' has to be > 0")

        if(is(neighbor, "ContNeighborhood")){
            bound <- 1 - 1/(2*E(Regressor, function(x, tau){pnorm(tau*abs(x))}, tau = risk@width))
            if(radius >= bound)
                stop("disjointness condition is violated!")
        }

        if(is(neighbor, "TotalVarNeighborhood")){
            bound <- E(Regressor, function(x, tau){pnorm(tau*abs(x))}, tau = risk@width) - 0.5
            if(radius >= bound)
                stop("disjointness condition is violated!")
        }

        c0 <- try(uniroot(getFixClipRegTS, lower = .Machine$double.eps^0.75, 
                    upper = upper, tol = tol, ErrorDistr = ErrorDistr, Regressor = Regressor, 
                    risk = risk, neighbor = neighbor)$root, silent = TRUE)

        Afct <- function(x, c0){
            if(x == 0) return(0)

            return(x^2*pnorm(c0/abs(x)))
        }
        A <- 1/(2*E(Regressor, Afct, c0 = c0) - E(Regressor, function(x){x^2}))

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        a <- -A*c0
        b <- 2*A*c0

        Risk <- getFiRiskRegTS(risk = risk, ErrorDistr = ErrorDistr, Regressor = Regressor, 
                            neighbor = neighbor, clip = c0, stand = A, sampleSize = sampleSize, 
                            Algo = Algo, cont = cont)

        return(list(A = as.matrix(A), a = a, b = b, d = NULL, risk = Risk, info = info))    
    })

setMethod("getFixRobRegTypeIC", signature(ErrorDistr = "Norm",
                                          Regressor = "UnivariateDistribution", 
                                          risk = "fiUnOvShoot", 
                                          neighbor = "CondNeighborhood"),
    function(ErrorDistr, Regressor, risk, neighbor, sampleSize, upper, maxiter, 
             tol, warn, Algo, cont){
        radiusCurve <- neighbor@radiusCurve
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

        radCx <- radiusCurve(x.vec)
        if(identical(all.equal(max(radCx), 0), TRUE))
            stop("'radiusCurve' has to be > 0")

        if(is(neighbor, "CondContNeighborhood")){
            test <- radCx - 1 + 1/(2*pnorm(risk@width*abs(x.vec)))
            if(!any(test < 0))
                stop("disjointness condition is violated!")
        }

        if(is(neighbor, "CondTotalVarNeighborhood")){
            test <- radCx - pnorm(risk@width*abs(x.vec)) + 0.5
            if(!any(test < 0))
                stop("disjointness condition is violated!")
        }

        b.vec <- numeric(length(x.vec))

        for(i in 1:length(x.vec)){
            if(test[i] >= 0)
                b.vec[i] <- 0
            else
                b.vec[i] <- try(uniroot(getFixClipRegTS, lower = .Machine$double.eps^0.75, 
                                upper = upper, tol = tol, ErrorDistr = ErrorDistr, 
                                Regressor = x.vec[i], risk = risk, neighbor = neighbor)$root, 
                            silent = TRUE)
        }

        if(is(Regressor, "DiscreteDistribution")){
            bfun <- function(x){ 
                ind <- (round(x, 8) == round(x.vec, 8))
                if(any(ind))
                    return(b.vec[ind])
                else
                    return(NA)
            }
        }else{
            if(is.finite(q(Regressor)(0)))
                yleft <- NA
            else
                yleft <- b[1]
            if(is.finite(q(Regressor)(1)))
                yright <- NA
            else
                yright <- b.vec[length(b.vec)]
            bfun <- approxfun(x.vec, b.vec, yleft = yleft, yright = yright)
        }

        Afct <- function(x, c0){ return(x^2*pnorm(c0(x))) }
        A <- 1/(2*E(Regressor, Afct, c0 = bfun) - E(Regressor, function(x){x^2}))

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getFiRiskRegTS(risk = risk, ErrorDistr = ErrorDistr, Regressor = Regressor, 
                            neighbor = neighbor, clip = bfun, stand = A, sampleSize = sampleSize, 
                            cont = cont)

        bfct <- function(x){bf <- bfun; A * bf(x)*abs(x)}
        body(bfct) <- substitute({bf <- bfun; A * bf(x)*abs(x)}, list(bfun = bfun, A = A))
        b <- RealRandVariable(Map = list(bfct), Domain = Reals())

        afct <- function(x){bf <- bfun; -A * bf(x)*abs(x)}
        body(afct) <- substitute({bf <- bfun; -A * bf(x)*abs(x)}, list(bfun = bfun, A = A))
        a <- RealRandVariable(Map = list(afct), Domain = Reals())

        return(list(A = as.matrix(A), a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
