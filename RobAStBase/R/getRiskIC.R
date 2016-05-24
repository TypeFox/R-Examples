###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "missing"),
    function(IC, risk, tol = .Machine$double.eps^0.25)
        getRiskIC(IC = IC, risk = risk,  L2Fam = eval(IC@CallL2Fam),
                  tol = tol))

setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, L2Fam, tol = .Machine$double.eps^0.25){
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(L2Fam@distribution)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        IC1 <- as(diag(dimension(IC@Curve)) %*% IC@Curve, "EuclRandVariable")

        bias <- E(L2Fam, IC1)
        Cov <- E(L2Fam, IC1 %*% t(IC1))

        prec <- checkIC(IC, L2Fam, out = FALSE)
        if(prec > tol)
            warning("The maximum deviation from the exact IC properties is ", prec,
                    "\nThis is larger than the specified 'tol' ",
                    "=> the result may be wrong")

        return(list(asCov = list(distribution = .getDistr(L2Fam), value = Cov - bias %*% t(bias))))
    })

###############################################################################
## trace of asymptotic covariance
###############################################################################
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "trAsCov",
                                 neighbor = "missing",
                                 L2Fam = "missing"),
    function(IC, risk, tol = .Machine$double.eps^0.25){
        getRiskIC(IC = IC, risk = risk,  L2Fam = eval(IC@CallL2Fam),
                  tol = tol)
    })

setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "trAsCov",
                                 neighbor = "missing",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, L2Fam, tol = .Machine$double.eps^0.25){
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(L2Fam@distribution)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        trCov <- getRiskIC(IC, risk = asCov(), L2Fam = L2Fam)$asCov
        trCov$value <- sum(diag(as.matrix(trCov$value)))

        prec <- checkIC(IC, L2Fam, out = FALSE)
        if(prec > tol)
            warning("The maximum deviation from the exact IC properties is ", prec,
                    "\nThis is larger than the specified 'tol' ",
                    "=> the result may be wrong")

        return(list(trAsCov = trCov))
    })

###############################################################################
## asymptotic Bias
###############################################################################
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "asBias",
                                 neighbor = "UncondNeighborhood",
                                 L2Fam = "missing"),
    function(IC, risk, neighbor, tol = .Machine$double.eps^0.25){
             getBiasIC(IC = IC, neighbor = neighbor, 
             biastype = biastype(risk), normtype = normtype(risk), tol = tol)
    })
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "asBias",
                                 neighbor = "UncondNeighborhood",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, neighbor, L2Fam, tol = .Machine$double.eps^0.25){
             getBiasIC(IC = IC, neighbor = neighbor, L2Fam = L2Fam, 
                       biastype = biastype(risk), normtype = normtype(risk), 
                       tol = tol)
    })
###############################################################################
## asymptotic MSE
###############################################################################
setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "asMSE",
                                 neighbor = "UncondNeighborhood",
                                 L2Fam = "missing"),
    function(IC, risk, neighbor, tol = .Machine$double.eps^0.25){
        L2Fam <- eval(IC@CallL2Fam)
        getRiskIC(IC = IC, risk = risk, neighbor = neighbor,
                  L2Fam = L2Fam, tol = tol)
    })

setMethod("getRiskIC", signature(IC = "IC",
                                 risk = "asMSE",
                                 neighbor = "UncondNeighborhood",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, neighbor, L2Fam, tol = .Machine$double.eps^0.25){
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(L2Fam@distribution)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        rad <- neighbor@radius
        if(rad == Inf) return(Inf)

        trCov <- getRiskIC(IC = IC, risk = trAsCov(), L2Fam = L2Fam)
        Bias <- getRiskIC(IC = IC, risk = asBias(), neighbor = neighbor, L2Fam = L2Fam)

        prec <- checkIC(IC, L2Fam, out = FALSE)
        if(prec > tol)
            warning("The maximum deviation from the exact IC properties is ", prec,
                    "\nThis is larger than the specified 'tol' ",
                    "=> the result may be wrong")

        nghb <- paste(neighbor@type, "with radius", neighbor@radius)

        return(list(asMSE = list(distribution = .getDistr(L2Fam),
                                 neighborhood = nghb,
                                 radius = neighbor@radius,
                                 value = trCov$trAsCov$value + rad^2*Bias$asBias$value^2)))
    })



