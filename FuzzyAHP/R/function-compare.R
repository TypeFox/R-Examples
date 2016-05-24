#' Function to mutually compare fuzzy data
#'
#' @description
#' This function compares fuzzy data \code{\linkS4class{FuzzyData}} to identify the optimal (best) cases.
#' The \code{\linkS4class{FuzzyData}} can only contain one fuzzy number, otherwise the comparison is not possible.
#' The calculation of type \code{"possibilityTheory"} can be time consuming as it is based on comparing each fuzzy
#' number to a maximal fuzzy number obtained from the set of fuzzy numbers.
#'
#'
#' @param fuzzyData A \code{\linkS4class{FuzzyData}}
#' @param type A \code{"character"} representing type of comparison. Currently implemented methods are
#' \code{"Chen"} and \code{"possibilityTheory"}.
#' @param progressBar logical value indicating if textual progress bar should be printed (default value \code{FALSE})
#'
#' @return A \code{"matrix"} of numeric value in case of \code{"Chen"} method or
#' a \code{"matrix"} with two values in case of \code{"possibilityTheory"}.
#'
#' @export
#' @rdname compareFuzzyNumbers-methods
#' @name compareFuzzyNumbers
setGeneric("compareFuzzyNumbers",
           signature = c("fuzzyData","type"),
           function(fuzzyData, type, progressBar = FALSE) standardGeneric("compareFuzzyNumbers"))

#' @rdname compareFuzzyNumbers-methods
#' @aliases compareFuzzyNumbers,FuzzyData,character-method
#' @importFrom utils setTxtProgressBar txtProgressBar
setMethod(
  f="compareFuzzyNumbers",
  signature(fuzzyData = "FuzzyData", type = "character"),
  definition=function(fuzzyData, type, progressBar)
  {

    if(ncol(fuzzyData@fnMin)>1){
      stop("Comparison of fuzzy values is only possible for datasets representing one fuzzy number.")
    }

    if(!(type == "Chen" || type == "possibilityTheory" )){
      stop("Type of comparion must be specified to either Chen or possibilityTheory!")
    }

    numberRows = length(fuzzyData@fnMin)

    if(type == "Chen"){
      minMin = min(fuzzyData@fnMin[,1], na.rm = TRUE)
      maxMax = max(fuzzyData@fnMax[,1], na.rm = TRUE)
      result = matrix(NA, nrow = numberRows, ncol = 1)
      result[,1] = (1/2)*(((fuzzyData@fnMax[,1]-minMin)/((maxMax - minMin)-(fuzzyData@fnModal[,1]-fuzzyData@fnMax[,1])))+1-((maxMax-fuzzyData@fnMin[,1])/((maxMax-minMin)+(fuzzyData@fnModal[,1]-fuzzyData@fnMin[,1]))))
    }
    else if(type == "possibilityTheory"){

      numberRows = nrow(fuzzyData@fnModal)
      result = matrix(NA, nrow = numberRows, ncol = 2)

      if(progressBar){
        pb = txtProgressBar(min=0, max=1, initial = 0, style = 3)
      }

      gmaxMin = max(fuzzyData@fnMin[,1], na.rm = TRUE)
      gmaxModal = max(fuzzyData@fnModal[,1], na.rm = TRUE)
      gmaxMax = max(fuzzyData@fnMax[,1], na.rm = TRUE)

      otherMax = FALSE

      if(length(which(fuzzyData@fnMin == gmaxMin)) < 2 &&
         length(which(fuzzyData@fnModal == gmaxModal)) < 2 &&
         length(which(fuzzyData@fnMax == gmaxMax)) < 2){
        otherMax = TRUE
      }

      for(i in 1:numberRows){

        if(is.na(fuzzyData@fnMin[i,]) || is.na(fuzzyData@fnModal[i,]) || is.na(fuzzyData@fnMax[i,])){
          result[i,1] = NA
          result[i,2] = NA
        }
        else{

          if(otherMax && (fuzzyData@fnMin[i,] == gmaxMin || fuzzyData@fnModal[i,] == gmaxModal || fuzzyData@fnMax[i,] == gmaxMax)){
            maxMin = max(fuzzyData@fnMin[-i,1], na.rm = TRUE)
            maxModal = max(fuzzyData@fnModal[-i,1], na.rm = TRUE)
            maxMax = max(fuzzyData@fnMax[-i,1], na.rm = TRUE)
          }
          else{
            maxMin = gmaxMin
            maxModal = gmaxModal
            maxMax = gmaxMax
          }

          # possibility
          if(fuzzyData@fnMax[i,] <= maxMin){
            result[i,1] = 0
          }
          else if(fuzzyData@fnModal[i,] >= maxModal){
            result[i,1] = 1
          }
          else{
            result[i,1] = .intersectionY(
              fuzzyData@fnMax[i,], 0,
              fuzzyData@fnModal[i,], 1,
              maxModal, 1,
              maxMin, 0)
          }

          # necessity
          if(fuzzyData@fnModal[i,] <= maxMin){
            result[i,2] = 0
          }
          else if(fuzzyData@fnMin[i,] >= maxModal){
            result[i,2] = 1
          }
          else{
            result[i,2] = .intersectionY(
              fuzzyData@fnModal[i,], 0,
              fuzzyData@fnMin[i,], 1,
              maxModal, 1,
              maxMin, 0
            )
          }

          if(progressBar){
            setTxtProgressBar(pb, (i/numberRows))
          }
        }
      }
      colnames(result) = c("possiblity", "necessity")
    }
    return(result)

  }
)





#' Function to rank results
#'
#' @description
#' This function ranks data from the highest value to the lowest. Essentially it does
#' \code{(nrow(data)+1) - rank(data, na.last ="keep", ties.method= "max")}.
#'
#'
#' @param data A matrix with one column. Usually an output of function \code{calculateAHP.}
#'
#' @return A \code{"matrix"} of numeric value that indicates ranking of each row, with the highest value ranked as 1.
#'
#' @export
#' @rdname compareResults-methods
#' @name compareResults
setGeneric("compareResults",
           function(data) standardGeneric("compareResults"))

#' @rdname compareResults-methods
#' @aliases compareResults,FuzzyData,character-method
#' @importFrom utils setTxtProgressBar txtProgressBar
setMethod(
  f="compareResults",
  signature(data = "matrix"),
  definition=function(data)
  {
    if(ncol(data)>1){
      stop("Comparing results is only possible for matrixes with one column.")
    }
    m = cbind((nrow(data)+1) - rank(data, na.last = "keep", ties.method= "max"))
    return(m)
  }
)





setGeneric(".intersectionY",
          function(x1,y1,x2,y2,x3,y3,x4,y4) standardGeneric(".intersectionY"))

setMethod(
  ".intersectionY",
  signature(x1 = "numeric", y1 = "numeric", x2 = "numeric", y2 = "numeric",
            x3 = "numeric", y3 = "numeric", x4 = "numeric", y4  = "numeric"),
  function (x1,y1,x2,y2,x3,y3,x4,y4)
  {
    x12 = x1 - x2
    x34 = x3 - x4
    y12 = y1 - y2
    y34 = y3 - y4

    c = x12 * y34 - y12 * x34

    a = x1 * y2 - y1 * x2
    b = x3 * y4 - y3 * x4

    #x = (a * x34 - b * x12) / c
    y = (a * y34 - b * y12) / c

    return(y)
  }
)
