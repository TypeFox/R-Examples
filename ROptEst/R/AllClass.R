.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
#    require("distr", character = TRUE, quietly = TRUE)
#    require("distrEx", character = TRUE, quietly = TRUE)
#    require("RandVar", character = TRUE, quietly = TRUE)
#    require("distrMod", character = TRUE, quietly = TRUE)
#    require("RobAStBase", character = TRUE, quietly = TRUE)
}


## asymptotic Anscombe risk
setClass("asAnscombe", representation(eff = "numeric"),
            prototype = prototype(eff = .95,
                             type = "optimal bias robust IC for given ARE in the ideal model"),
            contains = "asRiskwithBias",
            validity = function(object){
                if(any(object@eff <= 0|object@eff > 1))
                    stop("'eff' has to be in (0,1]")
                else TRUE
            })


## asymptotic L4 error
setClass("asL4", contains = "asGRisk",
            prototype = prototype(type = "asymptotic mean power 4 error"))
## asymptotic L1 error
setClass("asL1", contains = "asGRisk",
            prototype = prototype(type = "asymptotic mean absolute error"))
