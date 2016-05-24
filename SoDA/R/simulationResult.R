setClassUnion("randomGeneratorState", "integer")

setClass("simulationResult",
         representation( firstState = "randomGeneratorState",
                        expr = "language",
                        result = "ANY",
                        lastState = "randomGeneratorState"))

simulationResult <- function(value, seed) {
    x <- new("simulationResult")
    if(!missing(seed))
      set.seed(seed)
    else if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      stop("No initialized state for the random generator and no seed supplied")
    x@firstState <- .Random.seed
    x@expr <- substitute(value)
    x@result <- value
    x@lastState <- .Random.seed
    x
}

resetSeed <- function(object, last = FALSE) {
    if(!is(object, "simulationResult"))
      stop(
           gettextf("Supplied object must be a simulationResult, got class \"%s\"",
                    class(object)), domain = NA)
    assign(".Random.seed", (if(last) object@lastState else object@firstState),
           envir = .GlobalEnv
           )
}

setMethod("show", "randomGeneratorState",
          function(object) {
              if(is.numeric(object) && length(object) > 10) {
                  x <- format(object[-(1:2)])
                  cat("Uniform generator: ", object[[1]], "; Normal: ", object[[2]],
                      "\nSeeds: \n", sep = "")
                  cat(c(x[1:2], " ... ", x[length(x)]), sep = ", ")
                  cat("\n     (length =", length(x), ")\n")
              }
              else
                callNextMethod()
          })

setMethod("show", "simulationResult",
          function(object) {
              cat("Initial state: "); show(object@firstState)
              cat("\nFinal state:   "); show(object@lastState)
              cat("\nExpression: "); dput(object@expr)
              cat("\nSummary of result:\n"); print(summary(object@result))
          })
