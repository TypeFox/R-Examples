#' Generating recursive (acyclic) SEM models represented by
#' binary vectors. The generated models also include a model
#' from each model complexity.
#' @title Random SEM models.
#' @param nPop number of models to generate or population size.
#' @param numVar number of variables.
#' @param consMatrix \code{m by 2} binary \code{\link{matrix}}
#' representing constraint/prior knowledge,
#' where \code{m} is the number of constraint. For example, known that
#' variables 2 and 3 do not cause variable 1, then
#' \code{constraint <- matrix(c(2, 1, 3, 1), 2, 2, byrow=TRUE))} will be
#' the constraint matrix.
#' @param longitudinal \code{TRUE} for longitudinal data,
#' and \code{FALSE} for cross-sectional data.
#' @return \code{nPop} or \code{minPop} by \code{m} \code{\link{matrix}},
#' where \code{m} is the length of the binary vector depending
#' of the given number of variables
#' and also whether longitudinal or cross-sectional model.
#' @details This function generates \code{nPop} random SEM models which are
#' represented by binary vectors; 1 means there is a causal path from,
#' e.g., variable \code{A} to \code{B}
#' and 0 otherwise. In addition, the generated models
#' have passed the cyclic test to ensure they are all acyclic. The function
#' also includes \code{minPop} models which representing models
#' from each model complexity,
#' where \code{minPop = numVar(numVar-1)/2+1} if
#' \code{longitudinal = FALSE}, or
#' \code{minPop = (numVar(numVar-1)/2+1)+numVar^2} if
#' \code{longitudinal = TRUE}. If \code{nPop <= minPop} then
#' this function will generate \code{minPop} models.
#' @examples models <- modelPop(nPop=25, numVar=6,
#' longitudinal=FALSE, consMatrix = matrix(c(1, 2), 1, 2))
#' models
#' @author Ridho Rahmadi \email{r.rahmadi@cs.ru.nl}
#' @export
modelPop <- function(nPop=NULL, numVar = NULL,
                     longitudinal = NULL, consMatrix = NULL) {

  # check arguments
  if (!is.null(nPop)) {
    if (!is.numeric(nPop) || is.matrix(nPop)) {
      stop("Argument nPop should be positive numeric, e.g., 50.")
    }
  } else {
    nPop <- 50
  }

  if (!is.null(numVar)) {
    if (!is.numeric(numVar) || is.matrix(numVar)) {
      stop("Argument numVar should be positive numeric, e.g., 6.")
    }
  } else {
    stop("Argument numVar cannot be missing.")
  }

  if (!is.null(longitudinal)) {
    if (!is.logical(longitudinal)) {
      stop("Argument longitudinal should be either logical TRUE or FALSE.")
    }
  } else {
    stop("Argument longitudinal cannot be missing.")
  }

  if(!is.null(consMatrix)) {
    if(!is.matrix(consMatrix)) {
      stop("The constraints should be formed in a matrix.")
    }
  } else {
    stop("Argument consMatrix cannot be missing.")
  }

  #get consVector from consMatrix
  constraint1 <- convertCons(consMatrix, numVar)

  #get models from each model complexity
  if (longitudinal) {
    stringSize <- (numVar * numVar + (numVar * (numVar - 1)))
  } else {
    stringSize <- (numVar * (numVar - 1))
  }

  modelPopulation <- initialPopulation(numVar, stringSize,
                                       longitudinal, consMatrix)

  if (nPop <= nrow(modelPopulation)) {
    return(modelPopulation)
  } else {
    #complete the population
    modelPopulation <- rbind(modelPopulation,
                             genPopulation(nPop-nrow(modelPopulation), numVar,
                                           longitudinal, constraint1))

    return(modelPopulation)
  }

}
