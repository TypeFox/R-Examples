#' Compute the model \code{chi-square} and \code{model complexity}
#' of the given SEM models.
#' @title Scoring the given SEM models.
#' @param theData a data frame containing the data to which the model is to
#' be fit. If parameter \code{longitudinal} is \code{TRUE}, the data frame
#' should be reshaped such that the first \code{n} data points contain the
#' relations that occur in the first two time slices \code{t_0} and \code{t_1}.
#' The next \code{n} data points contain the relations that occur in
#' time slices \code{t_1} and \code{t_2}. The \code{i-th} subset of \code{n}
#' data points contain the relations in time slices \code{t_i-1} and \code{t_i}.
#' @param allModelString \code{m} by \code{n} \code{\link{matrix}} of
#' binary vectors representing models, where \code{m} is the number of models,
#' and \code{n} is the length of the binary vector.
#' @param numTime number of time slices. If a cross-sectional data then the
#' it is 1.
#' @param longitudinal \code{TRUE} for longitudinal data,
#' and \code{FALSE} for cross-sectional data.
#' @param co whether to use \code{"covariance"} or \code{"correlation"}
#' \code{\link{matrix}}.
#' @return a \code{\link{matrix}} of models including their fitness':
#' \code{chi-square} and \code{model complexity.}
#' @examples
#' \donttest{
#' the_data <- adhd
#' models <- modelPop(nPop=25, numVar=6, longitudinal=FALSE,
#' consMatrix = matrix(c(1, 2), 1, 2))
#'
#' model_fitness <- getModelFitness(theData=the_data,
#' allModelString=models, numTime=1, longitudinal=FALSE, co="covariance")
#' model_fitness}
#' @author Ridho Rahmadi \email{r.rahmadi@cs.ru.nl}
#' @export
getModelFitness <- function(theData = NULL, allModelString = NULL,
                            numTime = NULL, longitudinal = NULL, co = NULL) {

  # to check arguments
  # argument data
  if(!is.null(theData)) { # if data is supplied
    if (is.numeric(theData) && !(is.matrix(theData))) {
      stop("Data should be either a data frame or a matrix of numerical values.")
    } else if (!(is.numeric(theData)) && is.data.frame(theData)) {
      if (any(sapply(theData, is.numeric) == FALSE)) {
        stop("Data should be either a data frame or a matrix of numerical values.")
      }
    } else if (!is.numeric(theData)) {
      stop("Data should be either a data frame or a matrix of numerical values.")
    }
  } else { # if not supplied
    stop("Data cannot be missing")
  }

  if(!is.null(allModelString)) {
    if(!is.matrix(allModelString)) {
      stop("Argument allModelString should be formed in a matrix.")
    }
  } else {
    stop("Argument allModelString cannot be missing.")
  }

  if (!is.null(numTime)) {
    if (!is.numeric(numTime) || is.matrix(numTime)) {
      stop("Argument numTime should be positive numeric, e.g., 1 if cross-sectional data.")
    }
  } else {
    numTime <- 1
  }

  if (!is.null(longitudinal)) {
    if (!is.logical(longitudinal)) {
      stop("Argument longitudinal should be either logical TRUE or FALSE.")
    }
  } else {
    stop("Argument longitudinal cannot be missing.")
  }

  if (!longitudinal && numTime > 1) {
    stop("Cross-sectional data should have only one time slice, e.g., numTime = 1")
  } else if (longitudinal && numTime == 1) {
    stop("Longitudinal data should have more than one time slices, e.g., numTime = 2, with two time slices.")
  }

  if(!is.null(co)) {
    if(!is.character(co)) {
      stop("Argument co should be a string of characters, e.g., either covariance or correlation.")
    } else {
      covMatrix <- c("covariance", "correlation")
      if (!co %in% covMatrix) {
        stop("Argument co should be either covariance or correlation matrix.")
      }
    }
  } else {
    co <- "covariance"
  }


  #pre-steps
  if (longitudinal){
    #get the number of instances of one time slice
    numInstances <- nrow(theData) / (numTime - 1)

    #compute the size of Subset
    sizeSubset <- floor(numInstances / 2) * (numTime - 1)

    #compute sizeSubsetGetData, how many instances drawn from each time slice
    sizeSubsetGetData <- floor(numInstances / 2)

    #compute the number of variables,
    #asummed data already in reshaped network t_0..t_i
    numVar <- ncol(theData) / 2

    #compute the size of the a longitudinal model string
    stringSize <- (numVar * numVar + (numVar * (numVar - 1))) #inter + intra

    newData <- getDataLongi(theData, numTime, sizeSubsetGetData)

  } else {
    #for variable sizeSubset
    sizeSubset <- floor(nrow(theData) / 2)

    #for variable numVar, the number of variables in the data
    numVar <- ncol(theData)

    #the size of a model string
    stringSize <- (numVar * (numVar - 1))

    newData <- getDataCross(theData, sizeSubset)
  }


  #get the fitness'
  #imposing index l = 1
  allFitness <- gatherFitness(newData, allModelString, sizeSubset,
                                 numVar, 1, longitudinal, co)

  # remove the columns of index and BIC
  allFitness <- allFitness[, -c(3, 4)]

  return(cbind(allModelString, allFitness))

}
