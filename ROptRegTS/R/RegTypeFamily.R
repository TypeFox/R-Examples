## generating function
RegTypeFamily <- function(name, distribution = LMCondDistribution(), distrSymm, 
                          ErrorDistr = Norm(), ErrorSymm, main = 0, nuisance, trafo, 
                          param, props = character(0), RegDistr = Norm(), RegSymm,
                          Regressor = RealRandVariable(c(function(x){x}), Domain = Reals())){
    if(missing(name))
        name <- "regression type family"

    if(missing(distrSymm)) distrSymm <- NoSymmetry()

    if(missing(param)) 
        param <- ParamFamParameter(name = paste("parameter of", name), 
                            main = main, nuisance = nuisance, trafo = trafo)

    if(missing(ErrorSymm)) ErrorSymm <- NoSymmetry()

    if(missing(RegSymm)) RegSymm <- NoSymmetry()
 
    return(new("RegTypeFamily", name = name, distribution = distribution, 
               distrSymm = distrSymm, ErrorDistr = ErrorDistr, RegDistr = RegDistr, 
               RegSymm = RegSymm, Regressor = Regressor, param = param, props = props))
}

## access methods
setMethod("ErrorDistr", "RegTypeFamily", function(object) object@ErrorDistr)
setMethod("ErrorSymm", "RegTypeFamily", function(object) object@ErrorSymm)
setMethod("RegDistr", "RegTypeFamily", function(object) object@RegDistr)
setMethod("RegSymm", "RegTypeFamily", function(object) object@RegSymm)
setMethod("Regressor", "RegTypeFamily", function(object) object@Regressor)
