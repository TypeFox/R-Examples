###############################################################################
## get optimally robust IC for finite-sample under-/overshoot risk
###############################################################################
setMethod("getFixRobIC", signature(Distr = "Norm", 
                                   risk = "fiUnOvShoot", 
                                   neighbor = "UncondNeighborhood"),
    function(Distr, risk, neighbor, sampleSize, upper, maxiter, tol, warn, 
             Algo, cont){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            stop("'radius' has to be > 0")
        }
        
        if(is(neighbor, "ContNeighborhood"))
            if(radius >= 1 - 1/(2*pnorm(risk@width)))
                stop("disjointness condition is violated!")

        if(is(neighbor, "TotalVarNeighborhood"))
            if(radius >= pnorm(risk@width) - 0.5)
                stop("disjointness condition is violated!")
        
        c0 <- try(uniroot(getFixClip, lower = .Machine$double.eps^0.75, 
                    upper = upper, tol = tol, Distr = Distr, risk = risk, 
                    neighbor = neighbor)$root, silent = TRUE)
        A <- 1/(2*pnorm(c0)-1)

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        a <- -A*c0
        b <- 2*A*c0

        Risk <- getFiRisk(risk = risk, Distr = Distr, neighbor = neighbor, 
                          clip = c0, stand = A, sampleSize = sampleSize, 
                          Algo = Algo, cont = cont)

        return(list(A = as.matrix(A), a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
