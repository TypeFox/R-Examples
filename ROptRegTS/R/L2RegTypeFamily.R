## generating function
L2RegTypeFamily <- function(name, distribution = LMCondDistribution(), distrSymm,
                             main = 0, nuisance, trafo, param, props = character(0), 
                             L2deriv = EuclRandVarList(EuclRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                                        Domain = EuclideanSpace(dimension=2),
                                                                        dimension = 1)), 
                             ErrorDistr = Norm(), ErrorSymm, RegDistr = Norm(), RegSymm, 
                             Regressor = RealRandVariable(Map = list(function(x){x}), Domain = Reals()),
                             ErrorL2deriv = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                                             Domain = Reals())), 
                             ErrorL2derivSymm, ErrorL2derivDistr, ErrorL2derivDistrSymm, FisherInfo){
    if(missing(name)) 
        name <- "L2 differentiable regression type family"
    
    if(missing(distrSymm)) distrSymm <- NoSymmetry()

    if(missing(param)) 
        param <- ParamFamParameter(name = paste("Parameter of", name), 
                                main = main, nuisance = nuisance, trafo = trafo)

    if(missing(ErrorSymm)) ErrorSymm <- NoSymmetry()

    if(missing(RegSymm)) RegSymm <- NoSymmetry()

    if(missing(ErrorL2derivDistr))
        ErrorL2derivDistr <- imageDistr(RandVar = ErrorL2deriv, distr = ErrorDistr)

    if(missing(ErrorL2derivSymm)){
        nrvalues <- numberOfMaps(ErrorL2deriv)
        SymmList <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            SymmList[[i]] <- NonSymmetric()
        ErrorL2derivSymm <- new("FunSymmList", SymmList)
    }

    if(missing(ErrorL2derivDistrSymm)){
        nrvalues <- length(ErrorL2derivDistr)
        SymmList <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            SymmList[[i]] <- NoSymmetry()
        ErrorL2derivDistrSymm <- new("DistrSymmList", SymmList)
    }

    if(missing(FisherInfo)){
        if(is(distribution, "UnivariateCondDistribution")){
            dims <- length(param)
            L2deriv1 <- as(diag(dims) %*% L2deriv, "EuclRandVariable")
            L2.L2 <- L2deriv1 %*% t(L2deriv1)
            res <- numeric(length(L2.L2))
            for(i in 1:length(L2.L2)){
                fct <- function(x, cond, f1){ f1(cbind(cond,x)) }
                res[i] <- E(RegDistr, .condE, D1 = distribution, fct = fct, 
                            f1 = L2.L2@Map[[i]])                
            }            
            FisherInfo <- PosDefSymmMatrix(matrix(res, nrow = dims))
        }else{
            stop("not yet implemented")
        }
    }else{
        FisherInfo <- PosDefSymmMatrix(FisherInfo)
    }
    
    new("L2RegTypeFamily", name = name, distribution = distribution, distrSymm = distrSymm, 
        param = param, props = props, L2deriv = L2deriv, ErrorDistr = ErrorDistr, 
        RegDistr = RegDistr, RegSymm = RegSymm, Regressor = Regressor, 
        ErrorL2deriv = ErrorL2deriv, ErrorL2derivSymm = ErrorL2derivSymm, 
        ErrorL2derivDistr = ErrorL2derivDistr, ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
        FisherInfo = FisherInfo)
}

## access methods
setMethod("L2deriv", "L2RegTypeFamily", function(object) object@L2deriv)
setMethod("ErrorL2deriv", "L2RegTypeFamily", function(object) object@ErrorL2deriv)
setMethod("ErrorL2derivSymm", "L2RegTypeFamily", function(object) object@ErrorL2derivSymm)
setMethod("ErrorL2derivDistr", "L2RegTypeFamily", function(object) object@ErrorL2derivDistr)
setMethod("ErrorL2derivDistrSymm", "L2RegTypeFamily", function(object) object@ErrorL2derivDistrSymm)
setMethod("FisherInfo", "L2RegTypeFamily", function(object) object@FisherInfo)

## check centering of L2 derivative and Fisher Information
setMethod("checkL2deriv", "L2RegTypeFamily", 
    function(L2Fam, out = TRUE){ 
        dims <- length(L2Fam@param)
        L2 <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        cent <- E(L2Fam, L2)
        if(out)
            cat("precision of centering:\t", cent, "\n")

        consist <- E(L2Fam, L2 %*% t(L2)) - as(L2Fam@FisherInfo, "matrix")
        if(out){
            cat("precision of Fisher information:\n")
            print(consist)
        }
        res <- max(abs(cent), abs(consist))
        names(res) <- "maximum deviation"
        
        return(res)
    })
