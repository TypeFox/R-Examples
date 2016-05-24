################################################################################
################################################################################
## author Till Junge <till.junge@altermail.ch>                                ##
##                                                                            ##
## Copyright (c) UNIL (Universite de Lausanne)                                ##
## NCCR - LIVES (National Centre of Competence in Research "LIVES -           ##
## Overcoming vulnerability: life course perspectives",                       ##
## <http://www.lives-nccr.ch/>)                                               ##
##                                                                            ##
## spacom is free software: you can redistribute it and/or modify it under    ##
## the terms of the GNU General Public License as published by the Free       ##
## Software Foundation, either version 2 of the License or any later version. ##
##                                                                            ##
## spacom is distributed in the hope that it will be useful, but WITHOUT ANY  ##
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  ##
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      ##
## details, see <http://www.gnu.org/licenses/>.                               ##
################################################################################
################################################################################


## provides the MLSpawExactObject needed for the analysis and
## performs all consistency checks
MakeMLSpawExactObject <- function(individual.level.data,
                                  context.id,
                                  formula,
                                  precise.data,
                                  verbose,
                                  obj=NULL) {

  ## create an empty MLSpawExactObject to fill
  if (is.null(obj)) {
    obj <- new("MLSpawExactObject")
  }

  ## make sure individual.level.data is a data.frame
  obj@individual.level.data <- checkType(individual.level.data, c("data.frame"),
                                         "individual.level.data")

  ## make sure precise.data is a data.frame
  if (is.null(precise.data)) {
    precise.data <- data.frame(unique(individual.level.data[[context.id]]))
    names(precise.data) <- context.id
  }
  obj@precise.data <- checkType(precise.data, c("data.frame"),
                                "precise.data")

  ## make sure context.id is a name in contextual.data
  obj@context.id <-
    checkContextId(context.id = context.id,
                   individual.level.data.names =
                     names(obj@individual.level.data),
                   precise.data.names = names(obj@precise.data))

  ## check that the context correspond upper and lower level
  checkContexts(obj@precise.data[[obj@context.id]],
                obj@individual.level.data[[obj@context.id]])

  ## make sure the formula isn't a string anymore
  obj@formula <- checkFormula(formula)

  ## check verbose flag
  obj@verbose <- check.flag(verbose, "verbose")

  return(obj)
}

summary.MLSpawExactOutput <- function(object,...) {
  print(summary(object@lme))
  message("\nStandardised fixed effects:")
  print(object@beta)}

setGeneric("summary")
setMethod("summary",
          signature=c("MLSpawExactOutput"),
          definition=summary.MLSpawExactOutput)

setMethod("fixef",
          signature=c("MLSpawExactOutput"),
          definition=function(object, ...) {lme4::fixef(object@lme, ...)})
setMethod("ranef",
          signature=c("MLSpawExactOutput"),
          definition=function(object, ...) {lme4::ranef(object@lme, ...)})
setMethod("VarCorr", "MLSpawExactOutput",
          function(x, sigma, rdig) VarCorr(x@lme, sigma, rdig))
setMethod("AIC", "MLSpawExactOutput",
          function(object, ..., k=2) AIC(object@lme, ..., k=k))
setMethod("BIC", "MLSpawExactOutput",
          function(object, ..., k=2) BIC(object@lme, ..., k=k))

printSmlObject <- function(x, is.print=TRUE, ...) {
  print(x@lme, ...)
  cat("\nStandardised fixed effects:\n")
  print(x@beta, ...)
  invisible(x)
}
## print() is S3 generic ==> define an S3 method {or none; show() is good enough}:
print.MLSpawExactOutput <- printSmlObject

setMethod("show", signature="MLSpawExactOutput",
          function(object) printSmlObject(object, is.print=FALSE))

PerformMLSpawExact <- function(obj, ...) {
  ## create an output object
  output.obj <- new("MLSpawExactOutput")

  ## chop up the formula in order to manipulate it
  # # formula.str<- paste(as.character(obj@ formula[2]), '~',
  # #                     as.character(obj@ formula[3]))
  has.na.individuals <- any(is.na(obj@individual.level.data))
  has.na.contextuals <- any(is.na(obj@precise.data))
  if (has.na.contextuals || has.na.individuals) {
    warning(
      "There are NA's in your data!\n. This is almost certainly *not* what ",
      "you wanted. Spacom may not be able to compute standardized ",
      "coefficients correctly (although it will try).\n")
  }
  merged.data <- merge(obj@ individual.level.data,
                       obj@precise.data, by=obj@context.id)

  # #fixed.effect.formula = as.formula(formula.str)
  fixed.effect.formula <- as.formula(obj@formula)

  output.obj@lme <- lmer(formula=fixed.effect.formula,
                         data=merged.data, ...)

  ## compute the standardised coefficients for contextual data
  coefficient.names <- names(fixef(output.obj@lme))[-1]
  nb.coefficients <- length(coefficient.names)
  output.obj@beta <- numeric(nb.coefficients)
  names(output.obj@beta) <- coefficient.names
  cont.names <- names(obj@precise.data)
  sd.outcome <- sd(obj@individual.level.data[[as.character(obj@formula[2])]])
  for (name in coefficient.names) {
    if (name %in% cont.names) {
      sd.explanatory <- sd(as.numeric(obj@precise.data[[name]]))
    } else {
      sd.explanatory <- sd(as.numeric(obj@individual.level.data[[name]]))
    }
    output.obj@beta[[name]] <- fixef(output.obj@lme)[[name]]*sd.explanatory/sd.outcome
  }
  return(output.obj)
}

MLSpawExact <- function(individual.level.data,
                        context.id,
                        formula,
                        precise.data=NULL,
                        verbose=TRUE,
                        ...) {
  ## build the MLSpawExactobject
  obj <-
    MakeMLSpawExactObject(individual.level.data = individual.level.data,
                          context.id = context.id,
                          formula = formula,
                          precise.data = precise.data,
                          verbose=verbose)
  model <- PerformMLSpawExact(obj, ...)
  return(model)
}

