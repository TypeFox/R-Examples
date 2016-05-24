###############################################################################
## asymptotic Bias for various types
###############################################################################
setMethod("getBiasIC", signature(IC = "IC",
                                 neighbor = "UncondNeighborhood"),
    function(IC, neighbor, L2Fam, biastype = symmetricBias(),
             normtype = NormType(), tol = .Machine$double.eps^0.25,
             numbeval = 1e5){

        misF <- FALSE
        if(missing(L2Fam)){
            misF <- TRUE 
            L2Fam <- eval(IC@CallL2Fam)
        }
        D1 <- L2Fam@distribution
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(D1)))
            stop("dimension of 'Domain' of 'Curve' != dimension of 'img' of 'distribution' of 'L2Fam'")

        x <- as.matrix(r(D1)(numbeval))
        x <- as.matrix(x[!duplicated(x),])

        Bias <- .evalBiasIC(IC = IC, neighbor = neighbor, biastype = biastype,
                            normtype = normtype, x = x, trafo = trafo(L2Fam@param))

        prec <- if(misF) checkIC(IC, out = FALSE) else
                         checkIC(IC, L2Fam, out = FALSE)
        if(prec > tol)
            warning("The maximum deviation from the exact IC properties is ", prec,
                    "\nThis is larger than the specified 'tol' ",
                    "=> the result may be wrong")
        return(list(asBias = list(distribution = .getDistr(L2Fam),
                    neighborhood = neighbor@type, value = Bias)))
    })


### help functions ( not exported to namespace) for getRiskIC

setMethod(".evalBiasIC", signature(IC = "IC",
                                 neighbor = "ContNeighborhood",
                                 biastype = "BiasType"),
    function(IC, neighbor, biastype, normtype, x, trafo){
        ICx <- evalRandVar(as(diag(dimension(IC@Curve)) %*% IC@Curve,
                            "EuclRandVariable"),x)[,,1]
        return(max(fct(normtype)(ICx)))}
    )

setMethod(".evalBiasIC", signature(IC = "IC",
                                 neighbor = "TotalVarNeighborhood",
                                 biastype = "BiasType"),
    function(IC, neighbor, biastype, normtype, x, trafo){
        if(nrow(trafo) > 1)
            stop("not yet implemented for dimension > 1")
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        res <- evalRandVar(IC1, x)
        return(max(res) - min(res))}
    )

setMethod(".evalBiasIC", signature(IC = "IC",
                                 neighbor = "ContNeighborhood",
                                 biastype = "onesidedBias"),
    function(IC, neighbor, biastype, x, trafo){
        if(nrow(trafo) > 1)
            stop("not yet implemented for dimension > 1")
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        res <- evalRandVar(IC1, x)
        if (sign(biastype)>0)
             return(max(res))
        else return(-min(res))
    })

setMethod(".evalBiasIC", signature(IC = "IC",
                                 neighbor = "ContNeighborhood",
                                 biastype = "asymmetricBias"),
    function(IC, neighbor, biastype, x, trafo){
        if(nrow(trafo) > 1)
            stop("not yet implemented for dimension > 1")
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        res <- evalRandVar(IC1, x)
        return(max(res)/nu(biastype)[2] -
               min(res)/nu(biastype)[1])}
    )

