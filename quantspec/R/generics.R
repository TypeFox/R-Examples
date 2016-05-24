################################################################################
#' Generic functions for accessing attributes of objects
#'
#' These generic functions are needed to access the objects' attributes.
#' Note that the naming convention \code{getAttribute} was applied, where
#' \code{attribute} is the name of the attribute/slot of the class of the
#' object.
#'
#' @name generics-accessors
#'
#' @param object object from which to get the value
#' @param ... optional parameters; for documentation see the documentation of
#'             the methods to each of the generic.
#'
#' @seealso
#' For an overview on the classes of the framework, and all of their
#' attributes, see the class diagrams in the package description
#' [cf. \code{\link{quantspec-package}}].


## Class-FreqRep

#' @name generics-accessors
#' @aliases getY
#' @export
setGeneric("getY",
    function(object, ...){standardGeneric("getY")})

#' @name generics-accessors
#' @aliases getValues
#' @export
setGeneric("getValues",
    function(object, ...){standardGeneric("getValues")})

#' @name generics-accessors
#' @aliases getCoherency
#' @export
setGeneric("getCoherency",
    function(object, ...){standardGeneric("getCoherency")})

#' @name generics-accessors
#' @aliases getIsRankBased
#' @export
setGeneric("getIsRankBased",
    function(object, ...){standardGeneric("getIsRankBased")})

#' @name generics-accessors
#' @aliases getB
#' @export
setGeneric("getB",
    function(object, ...){standardGeneric("getB")})


## Class-LagEstimator

#' @name generics-accessors
#' @aliases getLagOperator
#' @export
setGeneric("getLagOperator",
    function(object, ...){standardGeneric("getLagOperator")})


## Class-LagOperator

#' @name generics-accessors
#' @aliases getMaxLag
#' @export
setGeneric("getMaxLag",
    function(object, ...){standardGeneric("getMaxLag")})


## Class-QRegEstimator

#' @name generics-accessors
#' @aliases getParallel
#' @export
setGeneric("getParallel",
    function(object, ...){standardGeneric("getParallel")})


## Class-QSpecQuantity

#' @name generics-accessors
#' @aliases getFrequencies
#' @export
setGeneric("getFrequencies",
    function(object, ...){standardGeneric("getFrequencies")})

#' @name generics-accessors
#' @aliases getLevels
#' @export
setGeneric("getLevels",
    function(object, ...){standardGeneric("getLevels")}) # j,


## Class-QuantileSD

#' @name generics-accessors
#' @aliases getMeanPG
#' @export
setGeneric("getMeanPG", function(object, ...){standardGeneric("getMeanPG")})

#' @name generics-accessors
#' @aliases getStdError
#' @export
setGeneric("getStdError", function(object, ...){standardGeneric("getStdError")})

#' @name generics-accessors
#' @aliases getN
#' @export
setGeneric("getN", function(object, ...){standardGeneric("getN")})

#' @name generics-accessors
#' @aliases getR
#' @export
setGeneric("getR", function(object, ...){standardGeneric("getR")})

#' @name generics-accessors
#' @aliases getType
#' @export
setGeneric("getType",
    function(object, ...){standardGeneric("getType")})

#' @name generics-accessors
#' @aliases getTs
#' @export
setGeneric("getTs", function(object, ...){standardGeneric("getTs")})


## Class-SmoothedPG

#' @name generics-accessors
#' @aliases getCoherencySdNaive
#' @export
setGeneric("getCoherencySdNaive", function(object, ...){standardGeneric("getCoherencySdNaive")})

#' @name generics-accessors
#' @aliases getSdNaive
#' @export
setGeneric("getSdNaive", function(object, ...){standardGeneric("getSdNaive")})

#' @name generics-accessors
#' @aliases getSdBoot
#' @export
setGeneric("getSdBoot", function(object, ...){standardGeneric("getSdBoot")})

#' @name generics-accessors
#' @aliases getPointwiseCIs
#' @export
setGeneric("getPointwiseCIs", function(object, ...){standardGeneric("getPointwiseCIs")})


## Class-Weight

#' @name generics-accessors
#' @aliases getDescr
#' @export
setGeneric("getDescr", function(object, ...){standardGeneric("getDescr")})


## Class-KernelWeight

#' @name generics-accessors
#' @aliases getW
#' @export
setGeneric("getW", function(object, ...){standardGeneric("getW")})

#' @name generics-accessors
#' @aliases getBw
#' @export
setGeneric("getBw", function(object, ...){standardGeneric("getBw")})

#' @name generics-accessors
#' @aliases getWnj
#' @export
setGeneric("getWnj", function(object, ...){standardGeneric("getWnj")}) #j,


################################################################################
#' Generic functions for implementation of methods of a class
#'
#' These generic functions need to be defined to allow for the automatic
#' dispaching mechanism.
#'
#' @name generics-functions
#'
#' @param object specifies the object from which the method is to be applied.
#' @param ... optional parameters; for documentation see the documentation of
#'             the methods to the generic.
#'
#' @seealso
#' For an overview on the classes of the framework, and all of their
#' methods, see the class diagrams in the package description
#' [cf. \code{\link{quantspec-package}}].

## Class-QuantileSD

#' @name generics-functions
#' @aliases increasePrecision
#' @export
setGeneric("increasePrecision",
    function(object, ...){standardGeneric("increasePrecision")})

## Class-BootPos

#' @name generics-functions
#' @aliases getPositions
#' @export
setGeneric("getPositions",
    function(object, ... ){standardGeneric("getPositions")}) # B=1

################################################################################
#' Generic functions for accessing associations of objects
#'
#' These generic functions are needed to access the objects' associated objects.
#' Note that the naming convention \code{getAssociatedObject} was applied, where
#' \code{AssociatedObject} is the name of the class of the associated object.
#'
#' @name generics-associations
#'
#' @param object object from which to get the associated object
#' @param ... optional parameters; for documentation see the documentation of
#'             the methods to each of the generic.
#'
#' @seealso
#' For an overview on the classes of the framework, and all
#' associations, see the class diagrams in the package description
#' [cf. \code{\link{quantspec-package}}].

## Class-QuantileSD
## Class-SmoothedPG

#' @name generics-associations
#' @aliases getQuantilePG
#' @export
setGeneric("getQuantilePG",
    function(object, ...){standardGeneric("getQuantilePG")})


## Class-FreqRep

#' @name generics-associations
#' @aliases getBootPos
#' @export
setGeneric("getBootPos",
    function(object, ...){standardGeneric("getBootPos")})

## Class-QuantilePG

#' @name generics-associations
#' @aliases getFreqRep
#' @export
setGeneric("getFreqRep",
    function(object, ...){standardGeneric("getFreqRep")})


## Class-IntegrQuantileSD

#' @name generics-associations
#' @aliases getQuantileSD
#' @export
setGeneric("getQuantileSD",
    function(object, ...){standardGeneric("getQuantileSD")})


## Class-SmoothedPG

#' @name generics-associations
#' @aliases getWeight
#' @export
setGeneric("getWeight",
    function(object, ...){standardGeneric("getWeight")})
