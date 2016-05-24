.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
### are the following calls to require necessary?
#    require("distr", character = TRUE, quietly = TRUE)
#    require("distrEx", character = TRUE, quietly = TRUE)
#    require("RandVar", character = TRUE, quietly = TRUE)
}

.onAttach <- function(library, pkg){
    unlockBinding(".distrModOptions", asNamespace("distrMod"))
    msga <- gettext(
    "Some functions from pkg's 'base' and 'stats' are intentionally masked ---see distrModMASK().\n"
                   )
    msgb <- gettext(
    "Note that global options are controlled by distrModoptions() ---c.f. ?\"distrModoptions\"."
                   )
    buildStartupMessage(pkg = "distrMod", msga, msgb,
                        library = library, packageHelp = TRUE,
        #                    MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
        VIGNETTE = gettext("There is a vignette to this package; try vignette(\"distrMod\").\n",
                           "Package \"distrDoc\" provides a vignette to the other distrXXX packages,",
                           "as well as to several related packages; try vignette(\"distr\").")
        )
    invisible()
}


distrModMASK <- function(library = NULL)
{
    infoShow(pkg = "distrMod", filename = "MASKING", library = library)
}


## norms to be used in prototype already
EuclideanNorm <- function(x){ 
    if(is.vector(x)) 
        return(abs(x))
    else 
        return(sqrt(colSums(x^2)))
}
QuadFormNorm <- function(x, A) sqrt(colSums(x*(A %*% x)))

################################
##
## Optional..-classes  Part I --- part II follows below
##
################################


## matrix or function or NULL -- a class for trafo's
setClassUnion("MatrixorFunction", c("matrix", "OptionalFunction"))
## matrix, numeric or NULL -- a class for covariance slots
setClassUnion("OptionalNumericOrMatrix", c("OptionalNumeric", "matrix"))
setClassUnion("OptionalNumericOrMatrixOrCall", c("OptionalNumericOrMatrix", "call"))
## DistrList or NULL -- a class for slot L2DerivDistr below
setClassUnion("OptionalDistrListOrCall", c("DistrList", "NULL", "call"))

################################
##
## classes
##
################################
### from Matthias' thesis / ROptEst

## symmetry of functions
setClass("FunctionSymmetry", contains = c("Symmetry", "VIRTUAL"))

## non-symmetric functions
setClass("NonSymmetric", contains = "FunctionSymmetry",
            prototype = prototype(type = "non-symmetric function",
                                  SymmCenter = NULL))

## even functions
setClass("EvenSymmetric", contains = "FunctionSymmetry",
            prototype = prototype(type = "even function",
                                  SymmCenter = numeric(0)))

## odd functions
setClass("OddSymmetric", contains = "FunctionSymmetry",
            prototype = prototype(type = "odd function",
                                  SymmCenter = numeric(0)))


## list of symmetry types
setClass(Class = "FunSymmList",
            prototype = prototype(list(new("NonSymmetric"))),
            contains = "list",
            validity = function(object){
                nrvalues <- length(object)
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "FunctionSymmetry"))
                        stop("element ", i, " is no 'FunctionSymmetry'")
                return(TRUE)
            })

### from Matthias' thesis / ROptEst
## Parameter of a parametric family of probability measures
setClass("ParamFamParameter",
            representation(main = "numeric",
                           nuisance = "OptionalNumeric",
                           fixed = "OptionalNumeric",
                           trafo = "MatrixorFunction"),
            prototype(name = "parameter of a parametric family of probability measures",
                      main = numeric(0), fixed = NULL, nuisance = NULL, 
                      trafo = new("matrix")),
            contains = "Parameter",
            validity = function(object){
                if(! ((is.matrix(object@trafo)) || is.function(object@trafo)))
                   stop("invalid transformation:\n",
                        "should be a matrix or a function")
                if(is.matrix(object@trafo)){
                ln.m <- length(object@main)
                ln.n <- length(object@nuisance)
                .validTrafo(object@trafo, ln.m, ln.m+ln.n) ### check validity
                return(TRUE)}
            })

setClass("ParamWithScaleFamParameter",
            contains = "ParamFamParameter")

setClass("ParamWithShapeFamParameter",
            representation(withPosRestr = "logical"),
            prototype(withPosRestr = TRUE),
            contains = "ParamFamParameter")

setClass("ParamWithScaleAndShapeFamParameter",
            contains = c("ParamWithScaleFamParameter",
                         "ParamWithShapeFamParameter")
         )

### from Matthias' thesis / ROptEst
## family of probability measures
setClass("ProbFamily", representation(name = "character",
                                      distribution = "Distribution",
                                      distrSymm = "DistributionSymmetry",
                                      props = "character"),
                       contains = "VIRTUAL")

### from Matthias' thesis / ROptEst
## parametric family of probability measures
setClass("ParamFamily",
            representation(param = "ParamFamParameter",
                           modifyParam = "function",
                           fam.call = "call",
                           startPar = "function",
                           makeOKPar = "function",
                           .withMDE = "logical",
                           .withEvalAsVar = "logical"
                           ### <- new !!! (not in thesis!)
                           ### a function with argument theta
                           ###  returning distribution P_theta
                           ),
            prototype(name = "parametric family of probability measures",
                      distribution = new("Norm"),
                      distrSymm = new("NoSymmetry"),
                      modifyParam = function(theta){ Norm(mean=theta) }, ### <- new !!! (not in thesis!)
                      fam.call = call("new", "ParamFamily"),
                      props = character(0),
                      makeOKPar = function(param)param,
                      startPar = function(x) {},
                      param = new("ParamFamParameter", main = 0, trafo = matrix(1)),
                      .withMDE = TRUE, .withEvalAsVar = TRUE),
            contains = "ProbFamily")


### from Matthias' thesis / ROptEst
## L2-differentiable parametric family of probability measures
setClass("L2ParamFamily",
            representation(L2deriv = "EuclRandVarList",
                           L2deriv.fct = "function", ## new: a function in theta which produces L2deriv
                           L2derivSymm = "FunSymmList",
                           L2derivDistr = "OptionalDistrListOrCall",
                           L2derivDistrSymm = "DistrSymmList",
                           FisherInfo = "PosSemDefSymmMatrix",
                           FisherInfo.fct = "function", ## new: a function in theta which produces FisherInfo
                           .withEvalL2derivDistr = "logical"
                           ),
            prototype(name = "L_2 differentiable parametric family of probability measures",
                      distribution = new("Norm"),
                      distrSymm = new("NoSymmetry"),
                      modifyParam = function(theta){ Norm(mean=theta) }, ### <- new !!! (not in thesis!)
                      fam.call = call("new", "L2ParamFamily"),
                      param = new("ParamFamParameter", main = 0, trafo = matrix(1)),
                      props = character(0),
                      startPar = function(x, ...) {},
                      L2deriv.fct = function(theta) {f <- function(x) {x-theta}
                                                     return(f)},
                      L2deriv = EuclRandVarList(RealRandVariable(Map = list(function(x)x),
                                                                 Domain = Reals())),
                      L2derivSymm = new("FunSymmList"),
                      L2derivDistr = UnivarDistrList(new("Norm")),
                      L2derivDistrSymm = new("DistrSymmList"),
                      FisherInfo.fct = function(theta)return(1),
                      FisherInfo = new("PosDefSymmMatrix"),
                      .withEvalL2derivDistr = TRUE),
            contains = "ParamFamily",
            validity = function(object){
                if(is(object@distribution, "UnivariateCondDistribution"))
                    stop("conditional distributions are not allowed in slot 'distribution'")

                if(!is(object@distrSymm, "NoSymmetry")){
                    if(!is(object@distrSymm@SymmCenter, "numeric"))
                        stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
                    if(length(object@distrSymm@SymmCenter) != dimension(img(object@distribution)))
                        stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
                }

                if(!is(object@FisherInfo, "PosDefSymmMatrix"))
                    warning("Fisher information is singular; not all parameter directions can be inferred.")

                dims <- length(object@param)
                if(ncol(object@FisherInfo) != dims)
                    stop(paste("dimension of 'FisherInfo' should be", dims))
                nrvalues <- numberOfMaps(object@L2deriv)
                if(nrvalues != length(object@L2derivSymm))
                    stop("number of Maps of 'L2deriv' != length of 'L2derivSymm'")
                if(nrvalues != length(object@L2derivDistr))
                    stop("number of Maps of 'L2deriv' != length of 'L2derivDistr'")
                if(nrvalues != length(object@L2derivDistrSymm))
                    stop("number of Maps of 'L2deriv' != length of 'L2derivDistrSymm'")
                if(dimension(Domain(object@L2deriv[[1]])) != dimension(img(object@distribution)))
                    stop("dimension of 'Domain' of 'L2deriv' != dimension of 'img' of 'distribution'")
                if(dimension(object@L2deriv) != dims)
                    stop("dimension of 'L2deriv' != dimension of parameters")

                return(TRUE)
            })

## L2-differentiable parametric family of probability measures generated by a group
setClass("L2GroupParamFamily",
            representation(LogDeriv = "function"),
            prototype(name = "L_2 differentiable parametric group family",
                      fam.call = call("new", "L2GroupParamFamily"),
                      Logderiv = function(x)x),
            contains = "L2ParamFamily")



## virtual in-between class for common parts in modifyModel - method
setClass("L2LocationScaleUnion",
            representation(locscalename = "character"),
         contains = c("L2GroupParamFamily","VIRTUAL")
        )

## virtual in-between class for common parts in modifyModel - method
setClass("L2ScaleShapeUnion", representation(scaleshapename ="character"),
         contains = c("L2GroupParamFamily","VIRTUAL")
        )

## virtual in-between class for common parts log/original scale methods
setClassUnion("L2ScaleUnion",
               c("L2LocationScaleUnion","L2ScaleShapeUnion")
               )

## L2-differentiable (univariate) location family
setClass("L2LocationFamily",
            prototype = prototype(locscalename = c("loc"="loc"),
                                  fam.call = call("new", "L2LocationFamily")),
            contains = "L2LocationScaleUnion")

## L2-differentiable (univariate) location scale family
setClass("L2ScaleFamily",
            prototype = prototype(locscalename = c("loc"="loc","scale"="scale"),
                                  fam.call = call("new", "L2ScaleFamily")),
            contains = "L2LocationScaleUnion")

## L2-differentiable (univariate) location and scale family
setClass("L2LocationScaleFamily",
            prototype = prototype(locscalename = c("loc"="loc","scale"="scale"),
                                  fam.call = call("new", "L2LocationScaleFamily")),
            contains = "L2LocationScaleUnion")


################################################################################
## Norm Classes
################################################################################

setClass("NormType", representation(name = "character", fct = "function"),
          prototype = prototype(name = "EuclideanNorm", fct = EuclideanNorm))

setClass("QFNorm", representation(QuadForm = "PosSemDefSymmMatrix"),
          prototype = prototype(name = "(Semi)Norm based on quadratic form",
                      QuadForm = new("PosSemDefSymmMatrix", matrix(1,1,1)),
                      fct = QuadFormNorm),
          contains = "NormType")

setClass("InfoNorm", contains = "QFNorm")
setClass("SelfNorm", contains = "QFNorm")

################################################################################
## Bias Classes
################################################################################

### from session 10-01-08 : class Bias type
setClass("BiasType", representation(name = "character"),
          contains = "VIRTUAL")

setClass("symmetricBias", prototype = prototype(name = "symmetric bias"),
          contains = "BiasType")

setClass("onesidedBias", representation(sign = "numeric"),
          prototype = prototype(name = "positive bias", sign = 1),
          contains = "BiasType")

setClass("asymmetricBias",
          representation(nu = "numeric"), ### weights acc. to paper
          prototype = prototype(name = "asymmetric bias", nu = c(1,1)),
          contains = "BiasType")

################################################################################
## Risk Classes
################################################################################

## risks (e.g., risk of estimator)
setClass("RiskType", representation(type = "character"),
          contains = "VIRTUAL")
## asymptotic risk
setClass("asRisk", contains = c("RiskType", "VIRTUAL"))
## asymptotic covariance
setClass("asCov", contains = "asRisk",
            prototype = prototype(type = "asymptotic covariance"))
## trace of asymptotic covariance
setClass("trAsCov", contains = "asRisk",
            prototype = prototype(type = "trace of asymptotic covariance"))

## asymptotic risk with bias

setClass("asRiskwithBias", representation(biastype = "BiasType",
                                          normtype = "NormType"),
          prototype = prototype(type = "asymptotic risk with bias",
                             biastype = new("symmetricBias"),
                             normtype = new("NormType", name = "EuclideanNorm",
                                             fct = EuclideanNorm)),
          contains = c("asRisk"))

## asymptotic Hampel risk
setClass("asHampel", representation(bound = "numeric"),
            prototype = prototype(bound = Inf,
                             type = "trace of asymptotic covariance for given bias bound"),
            contains = "asRiskwithBias",
            validity = function(object){
                if(any(object@bound <= 0))
                    stop("'bound' has to be positive")
                else TRUE
            })
## asymptotic bias
setClass("asBias",
            contains = "asRiskwithBias",
            prototype = prototype(type = "asymptotic bias"))

## convex asymptotic risk
setClass("asGRisk", contains = "asRiskwithBias")
## asymptotic mean square error
setClass("asMSE", contains = "asGRisk",
            prototype = prototype(type = "asymptotic mean square error"))
## asymptotic under-/overshoot probability
setClass("asUnOvShoot", representation(width = "numeric"),
            prototype = prototype(type = "asymptotic under-/overshoot probability"),
            contains = "asGRisk",
            validity = function(object){
                if(length(object@width) != 1)
                    stop("length of 'width' has to be 1")
                if(any(object@width <= 0))
                    stop("'width' has to be positive")
                else TRUE
            })
## finite-sample risk
setClass("fiRisk", contains = c("RiskType", "VIRTUAL"))
## finite-sample covariance
setClass("fiCov", contains = "fiRisk",
            prototype = prototype(type = "finite-sample covariance"))
## trace of finite-sample covariance
setClass("trFiCov", contains = "fiRisk",
            prototype = prototype(type = "trace of finite-sample covariance"))
## finite-sample Hampel risk
setClass("fiHampel", representation(bound = "numeric"),
            prototype = prototype(bound = Inf,
                             type = "finite-sample variance for given bias bound"),
            contains = "fiRisk",
            validity = function(object){
                if(any(object@bound <= 0))
                    stop("'bound' has to be positive")
                else TRUE
            })
## finite-sample mean square error
setClass("fiMSE", contains = "fiRisk",
            prototype = prototype(type = "finite-sample mean square error"))
## finite-sample bias
setClass("fiBias", contains = "fiRisk",
            prototype = prototype(type = "finite-sample bias"))
## finite-sample under-/overshoot probability
setClass("fiUnOvShoot", representation(width = "numeric"),
            prototype = prototype(type = "finite-sample under-/overshoot probability"),
            contains = "fiRisk",
            validity = function(object){
                if(length(object@width) != 1)
                    stop("length of 'width' has to be 1")
                if(any(object@width <= 0))
                    stop("'width' has to be positive")
                else TRUE
            })
## end Matthias' thesis
###############################################

setClass("asSemivar",
          contains = "asGRisk",
          prototype = prototype(type = "asymptotic Semivariance",
          biastype =  new("onesidedBias")))

#################################################
## "Estimate" classes
#################################################
setClass("Estimate",
         representation(name = "character",
                        estimate = "ANY",
                        samplesize = "numeric",
                        completecases = "logical",
                        asvar = "OptionalNumericOrMatrixOrCall",
                        Infos = "matrix",
                        estimate.call = "call",
                        nuis.idx = "OptionalNumeric",
                        fixed = "OptionalNumeric",
                        trafo = "list",
                        untransformed.estimate = "ANY",
                        untransformed.asvar = "OptionalNumericOrMatrixOrCall"),
         prototype(name = "Estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   estimate.call = call("{}"),
                   asvar = NULL,
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(0))},
                                mat = matrix(1)), ### necessary for comparison with unit matrix
                   nuis.idx = NULL,
                   fixed = NULL,
                   untransformed.estimate = NULL,
                   untransformed.asvar = NULL),
         validity = function(object){
            if(is.null(object@untransformed.estimate)){
               if(is.null(dim(object@estimate)))
                  len <- length(object@estimate)
               else
                  len <- dim(object@estimate)[1]
            }else{
               if(is.null(dim(object@untransformed.estimate)))
                  len <- length(object@untransformed.estimate)
               else
                  len <- dim(object@untransformed.estimate)[1]
            }
            if(!is.character(object@Infos))
                stop("'Infos' contains no matrix of characters")
            if(ncol(object@Infos)!=2)
                stop("'Infos' must have two columns")
            if(!is.null(object@nuis.idx))
                {if(any(object@nuis.idx<0) || any(object@nuis.idx>len))
                   stop(gettextf("'nuis.idx' must be in 1:%d", len))}
            else TRUE
         })

setClass("MCEstimate",
         representation(criterion = "numeric",
                        criterion.fct = "function",
                        method = "character",
                        optimwarn = "character",
                        startPar = "ANY"),
         prototype(name = "Minimum criterion estimate",
                   estimate = numeric(0),
                   samplesize = numeric(0),
                   completecases = logical(0),
                   asvar = NULL,
                   estimate.call = call("{}"),
                   criterion.fct =  function(){},
                   criterion = numeric(0),
                   method = "",
                   Infos = matrix(c(character(0),character(0)), ncol=2,
                                  dimnames=list(character(0), c("method", "message"))),
                   nuis.idx = NULL,
                   trafo = list(fct = function(x){
                                      list(fval = x, mat = matrix(1))},
                                mat = matrix(1)),
                   optimwarn = "",
                   startPar = NULL
                   ),
         contains = "Estimate")

## To Do: class MLEstimate which is compatible with class
## mle or maybe class summary.mle of package "stats4"

#################################################
## "Confint" classes
#################################################

setClass("Confint",
         representation(type = "character",
                        confint = "array",
                        call.estimate = "call",
                        name.estimate = "character",
                        samplesize.estimate = "numeric",
                        completecases.estimate = "logical",
                        trafo.estimate = "list",
                        nuisance.estimate = "OptionalNumeric",
                        fixed.estimate = "OptionalNumeric"
                        ),
         prototype(type = "",
                   confint = array(0),
                   call.estimate = call("{}"),
                   samplesize.estimate = numeric(0),
                   completecases.estimate = logical(0),
                   name.estimate = "",
                   trafo.estimate = list(fct = function(x){
                                             list(fval = x, mat = matrix(1))},
                                         mat = matrix(1)), ### necessary for comparison with unit matrix
                   nuisance.estimate = NULL,
                   fixed.estimate = NULL)
         )

################################
##
## Optional..-classes  Part II
##
################################

##  used for one common print method:
setClassUnion("ShowDetails", c("Estimate", "Confint",
                                "PosSemDefSymmMatrix", "ParamFamily",
                                "ParamFamParameter"))

