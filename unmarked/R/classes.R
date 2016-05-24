setClassUnion("optionalDataFrame", c("data.frame","NULL"))


setClassUnion("optionalMatrix", c("matrix","NULL"))


setClassUnion("optionalNumeric", c("numeric","NULL"))


setClassUnion("optionalCharacter", c("character","NULL"))

setClassUnion("optionalList", c("list","NULL"))

# Compute standard error of an object.
setGeneric("SE",
    def = function(obj, ...) {
      standardGeneric("SE")
    })

setGeneric("plot")
setGeneric("predict")
setGeneric("vcov")
setGeneric("coef")
setGeneric("summary")
setGeneric("update")
setGeneric("confint")
setGeneric("profile")
setGeneric("head")
setGeneric("fitted")
setGeneric("simulate")
setGeneric("residuals")
setGeneric("hist")
setGeneric("logLik")


# Compute linear combinations of parameters.
setGeneric("linearComb",
    function(obj, coefficients, ...) {
      standardGeneric("linearComb")
    })

# Transform an object to it's natural scale.
setGeneric("backTransform",
    function(obj, ...) {
      standardGeneric("backTransform")
    })


setGeneric("hessian",	function(object) standardGeneric("hessian"))


setGeneric("LRT", function(m1, m2) standardGeneric("LRT"))


# TODO: make parent class unmarkedEstimate
# TODO: make binomial detection child class
# TODO: make binomial occ child class
# TODO: make poisson abundance child class
# TODO: make show method for each of these classes.
# TODO: make unmarkedFit show that calls the respective children.
# TODO: separate unmarkedFit class for each model type?... would contain different estimate types
