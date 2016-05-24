.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
#    require("distr", character = TRUE, quietly = TRUE)
#    require("distrEx", character = TRUE, quietly = TRUE)
#    require("distrMod", character = TRUE, quietly = TRUE)
#    require("RandVar", character = TRUE, quietly = TRUE)
}

.onAttach <- function(library, pkg){
    unlockBinding(".RobAStBaseOptions", asNamespace("RobAStBase"))
    msga <- gettext(
    "Some functions from pkg's 'stats' and 'graphics' are intentionally masked ---see RobAStBaseMASK().\n"
                   )
    msgb <- gettext(
    "Note that global options are controlled by RobAStBaseoptions() ---c.f. ?\"RobAStBaseoptions\"."
                   )
    buildStartupMessage(pkg = "RobAStBase", msga, msgb,
                        library = library, packageHelp = TRUE
        #                    , MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf"
        #                    , VIGNETTE = gettext("Package \"distrDoc\" provides a vignette to this package as well as to several related packages; try vignette(\"distr\").")
        )
    invisible()
}

RobAStBaseMASK <- function(library = NULL)
{
    infoShow(pkg = "RobAStBase", filename = "MASKING", library = library)
}

## neighborhood
setClass("Neighborhood",
            representation(type = "character",
                           radius = "numeric"), 
            contains = "VIRTUAL")
## unconditional (errors-in-variables) neighborhood
setClass("UncondNeighborhood", contains = c("Neighborhood", "VIRTUAL"))
## unconditional convex contamination neighborhood
setClass("ContNeighborhood", contains = "UncondNeighborhood",
            prototype = prototype(type = "(uncond.) convex contamination neighborhood",
                                  radius = 0))
## unconditional total variation neighborhood
setClass("TotalVarNeighborhood", contains = "UncondNeighborhood",
            prototype = prototype(type = "(uncond.) total variation neighborhood",
                                  radius = 0))
## robust model
setClass("RobModel",
            representation(center = "ProbFamily",
                           neighbor = "Neighborhood"),
            contains = "VIRTUAL")
## robust model with fixed (unconditional) neighborhood
setClass("FixRobModel",
            prototype = prototype(center = new("ParamFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel",
            validity = function(object){
                if(!is(object@neighbor, "UncondNeighborhood"))
                    stop("'neighbor' is no unconditional neighborhood")
                if(any(object@neighbor@radius < 0 || object@neighbor@radius > 1))
                    stop("neighborhood radius has to be in [0, 1]")
                else return(TRUE)
            })
## robust model with infinitesimal (unconditional) neighborhood
setClass("InfRobModel",
            prototype = prototype(center = new("L2ParamFamily"),
                                  neighbor = new("ContNeighborhood")),
            contains = "RobModel",
            validity = function(object){
                if(!is(object@neighbor, "UncondNeighborhood"))
                    stop("'neighbor' is no unconditional neighborhood")
                if(any(object@neighbor@radius < 0))
                    stop("'radius' has to be in [0, Inf]")
                else return(TRUE)
            })
## Weights
setClass("RobAStControl", representation(name ="character"),
          contains = "VIRTUAL")

setClass("RobWeight", representation(name = "character", weight = "function"), 
          prototype(name = "some weight", weight = function(x) 1))
setClass("BoundedWeight", representation(clip = "numeric"), 
          prototype(clip = 1), contains = "RobWeight")
setClass("BdStWeight", representation(stand = "matrix"), 
          prototype(stand = matrix(1)), contains = "BoundedWeight")
setClass("HampelWeight", representation(cent = "numeric"), 
          prototype(cent = 0), contains = "BdStWeight")




## Influence curve/function with domain: EuclideanSpace
setClass("InfluenceCurve", 
            representation(name = "character", 
                           Curve = "EuclRandVarList", 
                           Risks = "list",
                           Infos = "matrix"),
            validity = function(object){
                if(!is(Domain(object@Curve[[1]]), "EuclideanSpace"))
                    stop("The domain of 'Curve' has to be a Euclidean space")
                if(!is.character(object@Infos))
                    stop("'Infos' contains no matrix of characters")
                for(char in names(object@Risks))
                    if(!extends(char, "RiskType"))
                        stop(paste(char, "is no valid 'RiskType'"))
                if(ncol(object@Infos)!=2)
                    stop("'Infos' must have two columns")
                else TRUE
            })
## partial incluence curve
setClass("IC", representation(CallL2Fam = "call",
                              modifyIC = "OptionalFunction"),
            prototype(name = "square integrable (partial) influence curve",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                               Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      modifyIC = NULL),
            contains = "InfluenceCurve",
            validity = function(object){
                L2Fam <- eval(object@CallL2Fam)
                trafo <- trafo(L2Fam@param)
                if(nrow(trafo) != dimension(object@Curve))
                    stop("wrong dimension of 'Curve'")
                if(dimension(Domain(L2Fam@L2deriv[[1]])) != dimension(Domain(object@Curve[[1]])))
                    stop("dimension of 'Domain' of 'L2deriv' != dimension of 'Domain' of 'Curve'")

                return(TRUE)
            })
## HampIC -- common mother class to ContIC and TotalVarIC 
setClass("HampIC", 
            representation(stand = "matrix",
                           lowerCase = "OptionalNumeric",
                           neighborRadius = "numeric",
                           weight = "RobWeight",
                           biastype = "BiasType",
                           normtype = "NormType"), 
            prototype(name = "IC of total-var or contamination type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                    Domain = Reals())),
                      Risks = list(),  weight = new("RobWeight"),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      modifyIC = NULL,
                      stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0, 
                      biastype = symmetricBias(), 
                      NormType = NormType()),
            contains = "IC",
            validity = function(object){
                if(any(object@neighborRadius < 0)) # radius vector?!
                    stop("'neighborRadius' has to be in [0, Inf]")
                if(!is.null(object@lowerCase))
                    if(length(object@lowerCase) != nrow(object@stand))
                        stop("length of 'lowerCase' != nrow of standardizing matrix")
                L2Fam <- eval(object@CallL2Fam)
                if(!identical(dim(trafo(L2Fam@param)), dim(object@stand)))
                    stop(paste("dimension of 'trafo' of 'param' != dimension of 'stand'"))
                return(TRUE)
            })
## (partial) influence curve of contamination type
setClass("ContIC", 
            representation(clip = "numeric",
                           cent = "numeric"), 
            prototype(name = "IC of contamination type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                                    Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      modifyIC = NULL,
                      clip = Inf, cent = 0, stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0, weight = new("HampelWeight"),
                      biastype = symmetricBias(), NormType = NormType()),
            contains = "HampIC",
            validity = function(object){
                if(length(object@cent) != nrow(object@stand))
                    stop("length of centering constant != nrow of standardizing matrix")
                if((length(object@clip) != 1) && (length(object@clip) != length(object@Curve)))
                    stop("length of clipping bound != 1 and != length of 'Curve'")
                if(!is(weight,"HampelWeight")) 
                    stop("Weight has to be of class 'HampelWeight'")
                return(TRUE)
            })
## (partial) influence curve of total variation type
setClass("TotalVarIC",
            representation(clipLo = "numeric",
                           clipUp = "numeric"),
            prototype(name = "IC of total variation type",
                      Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}),
                                                               Domain = Reals())),
                      Risks = list(),
                      Infos = matrix(c(character(0),character(0)), ncol=2,
                                dimnames=list(character(0), c("method", "message"))),
                      CallL2Fam = call("L2ParamFamily"),
                      modifyIC = NULL,
                      clipLo = -Inf, clipUp = Inf, stand = as.matrix(1),
                      lowerCase = NULL,
                      neighborRadius = 0, weight = new("BdStWeight"),
                      biastype = symmetricBias(), NormType = NormType()),
            contains = "HampIC",
            validity = function(object){
                if((length(object@clipLo) != 1) && (length(object@clipLo) != length(object@Curve)))
                    stop("length of lower clipping bound != 1 and != length of 'Curve'")
                if((length(object@clipLo) != 1) && (length(object@clipLo) != length(object@Curve)))
                    stop("length of upper clipping bound != 1 and != length of 'Curve'")
                if(!is(weight,"BdStWeight")) 
                    stop("Weight has to be of class 'BdStWeight'")
                return(TRUE)
            })

## ALEstimate
setClassUnion("OptionalInfluenceCurve", c("InfluenceCurve", "NULL"))
setClassUnion("StartClass", c("numeric", "matrix", "function", "Estimate"))
setClass("pICList",
          prototype = prototype(list()),
            contains = "list",
            validity = function(object){
                nrvalues <- length(object)
                if(nrvalues){
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "OptionalInfluenceCurve"))
                        stop("element ", i, " is no 'OptionalInfluenceCurve'")
                }
                return(TRUE)
            })
setClassUnion("OptionalpICList", c("pICList", "NULL"))
setClass("ALEstimate",
         representation(pIC = "OptionalInfluenceCurve",
                        asbias = "OptionalNumeric"),
         prototype(name = "Asymptotically linear estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   estimate.call = call("{}"),
                   asvar = NULL,
                   asbias = NULL,
                   pIC = NULL,
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1)), ### necessary for comparison with unit matrix
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   completecases = logical(0),
                   untransformed.estimate = NULL,
                   untransformed.asvar = NULL),
         contains = "Estimate")
setClass("kStepEstimate", 
         representation(steps = "integer",
                        pICList = "OptionalpICList",
                        ICList = "OptionalpICList",
                        start = "StartClass",
                        startval = "matrix",
                        ustartval = "matrix",
                        ksteps = "OptionalMatrix",
                        uksteps = "OptionalMatrix"),
         prototype(name = "Asymptotically linear estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   estimate.call = call("{}"),
                   steps = integer(0),
                   asvar = NULL,
                   asbias = NULL,
                   pIC = NULL,
                   pICList = NULL,
                   ICList = NULL,
                   ksteps = NULL,
                   uksteps = NULL,
                   start = matrix(0),
                   startval = matrix(0),
                   ustartval = matrix(0),
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1)), ### necessary for comparison with unit matrix
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   untransformed.estimate = NULL,
                   untransformed.asvar = NULL),
         contains = "ALEstimate")
setClass("MEstimate", 
         representation(Mroot = "numeric"),
         prototype(name = "Asymptotically linear estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   estimate.call = call("{}"),
                   Mroot = numeric(0),
                   asvar = NULL,
                   asbias = NULL,
                   pIC = NULL,
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1)), ### necessary for comparison with unit matrix
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   untransformed.estimate = NULL,
                   untransformed.asvar = NULL),
         contains = "ALEstimate")
#################################################
## "cutoff" class
#################################################
setClass("cutoff", representation = representation(name = "character",
                                                   fct = "function",
                                                   cutoff.quantile = "numeric"),
                   prototype = prototype(name = "empirical",
                                         fct = function(data) quantile(data),
                                         cutoff.quantile = 0.95))


#################################################
# new risk classes
#################################################
setClass("interpolRisk", representation = representation(samplesize="numeric"),
                         contains = c("VIRTUAL", "RiskType"))
setClass("OMSRRisk", contains = "interpolRisk", prototype=prototype(type=".OMSE", samplesize=100))
setClass("RMXRRisk", contains = "interpolRisk", prototype=prototype(type=".RMXE", samplesize=100))
setClass("MBRRisk", contains = "interpolRisk", prototype=prototype(type=".MBRE",samplesize=100))
