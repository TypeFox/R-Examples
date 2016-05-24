###############################################################################
## Location M estimator (univariate location)
###############################################################################
setMethod("locMEstimator", signature(x = "numeric", IC = "InfluenceCurve"),
    function(x, IC, eps = .Machine$double.eps^0.5, na.rm = TRUE){
        es.call <- match.call()
        es.call[[1]] <- as.name("locMEstimator")
        if(numberOfMaps(IC@Curve) > 1)
            stop("number of Maps of 'IC' has to be 1")

        completecases <- complete.cases(x)
        if(na.rm) x <- na.omit(x)
        
        mest <- function(theta, x, IC){
            return(rowSums(evalIC(IC, as.matrix(x-theta))))
        }
        res <- uniroot(f = mest, interval = c(min(x), max(x)), 
             tol = eps, x = x, IC = IC)

        if(is(IC, "IC")){
            L2Fam <- eval(CallL2Fam(IC))
            Infos <- matrix(c("locMEstimator", 
                            paste("Location M estimate for", name(L2Fam))),
                            ncol = 2)
            colnames(Infos) <- c("method", "message")
            asVar <- getRiskIC(IC, risk = asCov(), L2Fam = L2Fam)$asCov$value
            asBias <- getRiskIC(IC, risk = asBias(), 
                                neighbor = ContNeighborhood(1), 
                                L2Fam = L2Fam)$asBias$value
                                 
            names(res$root) <- nms <- locscalename(L2Fam)
            asVar <- PosDefSymmMatrix(asVar)
            dimnames(asVar) <- list(nms, nms)
            names(asBias) <- nms
        }else{
            Infos <- matrix(c("locMEstimator", "Location M estimate"), ncol = 2)
            colnames(Infos) <- c("method", "message")
            asvar <- NULL
        }
        new("MEstimate", name = "Location M estimate", estimate = res$root, 
            completecases = completecases,
            estimate.call = es.call, pIC = IC, Mroot = res$f.root, 
            Infos = Infos, samplesize = length(x), asvar = asVar, asbias = asBias)
    })
