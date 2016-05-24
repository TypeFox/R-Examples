#' Repairing a SEM model that is cyclic.
#' @param stringModel binary vector with length
#' \code{n^2+(n(n-1))} if \code{longitudinal = TRUE},
#' or \code{n(n-1)} binary vector if
#' \code{FALSE}, where \code{n} is the number of variables (equal to argument
#' \code{numVar}).
#' @param numVar number of variables.
#' @param longitudinal \code{TRUE} for longitudinal data,
#' and \code{FALSE} for cross-sectional data.
#' @return a binary vector with the same length of input, representing a
#' repaired or acyclic model.
#' @examples
#' num_vars <- 6
#' longi_a <- FALSE
#' longi_b <- TRUE
#'
#' # Assume that the generated model below is cyclic
#' # a cross-sectional model
#' model_a <- round(runif(num_vars * num_vars))
#' # a longitudinal model
#' model_b <- c(round(runif(num_vars * num_vars)),
#' round(runif(num_vars * (num_vars-1))))
#'
#' repaired_model_a <- repairCyclicModel(stringModel=model_a, numVar=num_vars,
#' longitudinal=longi_a)
#'
#' repaired_model_b <- repairCyclicModel(stringModel=model_b, numVar=num_vars,
#' longitudinal=longi_b)
#'
#' repaired_model_a
#' repaired_model_b
#'
#' @details The main idea of this function is to seek cyclic(s) with
#' any possible length from a given model, and then to cut the cyclic,
#' so as to make the model acyclic. Moreover, this function is used in
#' \code{\link{stableSpec}} to ensure no cyclic model in the computation.
#' @author Ridho Rahmadi \email{r.rahmadi@cs.ru.nl}
#' @export
repairCyclicModel <- function(stringModel = NULL, numVar = NULL,
                              longitudinal = NULL) {

  if (!is.null(stringModel)) {
    if (!(all(stringModel %in% 0:1)) || is.matrix(stringModel)) {
      stop("Argument numVar should be binary vector.")
    }

  } else {
    stop("Argument stringModel cannot be missing.")
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

  # convert model string into matrix
  # if longitudinal = TRUE, then only the intraString would be checked
  theModel <- stringToMatrix1(stringModel, numVar, longitudinal)

  # if the model is cyclic
  if(!ggm::isAcyclic(theModel)) { # instead if(ggm::isAcyclic(theModel))

    #repair the model
    theModel <- cycleRepair(theModel)

    #to convert back from matrix to vector
    diag(theModel) <- NA
    stringModel_intra <- as.vector(theModel)

    if (longitudinal) {
      #after being repaired, the intraString is combined
      #again with its interString
      resModel <- c(stringModel[1:(numVar * numVar)],
                    stringModel_intra[!is.na(stringModel_intra)])
    } else { #if cross-sectional
      resModel <- stringModel_intra[!is.na(stringModel_intra)]
    }
    return(resModel)

  } else {
    # if Acyclic
    # both longitudinal and cross-sectional
    return(stringModel)
  }
}

cycleRepair <- function(theModel) {
  for (i in 1:nrow(theModel)) {

    squaredModel <- matrixcalc::matrix.power(theModel, i)
    cycle <- sum(diag(squaredModel))

    while (cycle != 0) {
      theIndex <- which(diag(squaredModel) != 0)

      #take the first and last index of theIndex
      first <- theIndex[1]
      last <- tail(theIndex, n=1)

      arc <- FALSE
      while(arc == FALSE) {
        k <- sample(first:last, 1)
        l <- sample(first:last, 1)

        if(theModel[k,l] != 0) {
          theModel[k,l] <- 0
          arc <- TRUE
        }
      }
      squaredModel <- matrixcalc::matrix.power(theModel, i)
      cycle <- sum(diag(squaredModel))
    }
  }
  return(theModel)
}
