setClass("target",
         representation(dimension="numeric", parameters="list",
                        type = "character", name = "character",
                        rinit = "function", dinit = "function",
                        generate = "function",
                        logdensity = "function", logdensityupdate =
                        "function", updateavailable = "logical"))
setGeneric("target", function(...)standardGeneric("target"))
target.constructor <- function(..., dimension, parameters, name,
                               type, rinit, dinit, logdensity, generate,
                               logdensityupdate){
  if (missing(dimension))
    stop(sQuote("dimension"), "has to be specified.")
  if (missing(parameters))
    parameters <- list()
  if (missing(name))
    name <- "unspecified"
  if (missing(type))
    type <- "continuous"
  if (missing(rinit))
    stop(sQuote("rinit"), "has to be specified to draw initial points of the Markov chains")
  if (missing(dinit)){
    cat("dinit has to be specified if you want to use SMC")
    dinit <- function(x) 0
  }
  if (missing(logdensity))
    stop(sQuote("logdensity"), "has to be specified (up to an (additive) normalizing constant)")    
  if (missing(generate))
      generate <- function(size, parameters, ...) stop(sQuote("generate"), "is not specified")
  if (missing(logdensityupdate)){
      logdensityupdate <- function(...) NULL
      updateavailable <- FALSE
  } else { updateavailable <- TRUE }

  new("target", dimension = dimension, parameters = parameters,
      name = name, type = type, rinit = rinit, dinit = dinit, logdensity = logdensity,
      generate = generate, logdensityupdate = logdensityupdate,
      updateavailable = updateavailable)
}
setMethod("target",
          definition = function(..., dimension, parameters, name, type,
                                rinit, dinit, logdensity, generate, logdensityupdate){
            target.constructor(dimension = dimension,
                               parameters = parameters, name = name,
                               type = type, rinit = rinit, dinit = dinit,
                               logdensity = logdensity,
                               generate = generate, logdensityupdate =
                               logdensityupdate)
          })

setMethod(f = "show", signature = "target",
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("*name:", object@name, "\n")
            cat("*type:", object@type, "\n")
            cat("*dimension of the state space:", object@dimension, "\n")
            cat("*target parameters:", names(object@parameters), "\n")
            cat("*update available:", object@updateavailable, "\n")
          })



