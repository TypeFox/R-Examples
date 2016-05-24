#' @include 00Utils.R
NULL

setClass("DataClass", representation(x="data.frame", y="factor", variance="logical", finite="logical", completeobs="logical", classbalance="logical", ntopratiotwoplus="logical", mindimensions="logical"))

validatedata <- function(object){

  data <- object

  temp <- length(caret::nearZeroVar(data@x))
  data@variance <- temp==0

  temp1 <- all(apply(data@x, 1:2, is.finite))
  data@finite <- temp1==TRUE

  temp2 <- any(apply(data@x, 1:2, is.na))
  data@completeobs <- temp2==FALSE

  temp3 <- length(caret::nearZeroVar(data.frame(data@y)))
  data@classbalance <- temp3==0

  temp5 <- nrow(data@x) > (2*ncol(data@x))
  data@ntopratiotwoplus <- temp5==TRUE

  temp6 <- all(dim(data@x) > c(20,3))
  data@mindimensions <- temp6

  minimumycheck <- min(table(data@y))
  if (minimumycheck < 5) {stop("One level in the factor variable has less than five observations.")}

  return(data)

}

#' initializedataclassobject
#'
#' initializedataclassobject is a function to construct a DataClass object.
#'
#' @param data (data.frame)
#' @details Argument 'data' must have only numeric columns and one factor column.
#' @export
#' @examples
#' ## dataobject <- initializedataclassobject(data.frame(x=scale(dataobject@@x), dataobject@@y))

initializedataclassobject <- function(data){

  if(class(data)!="data.frame"){stop("Argument 'data' must be a data frame.")}
  if(sum(sapply(data, is.factor)==TRUE)!=1) {stop("Argument 'data' must have one and only one factor column.")}
  if(sum(sapply(data, is.numeric)==TRUE)!=ncol(data)-1) {stop("Argument 'data' must have only numeric columns and one factor column.")}

  dataclassobject <- new("DataClass")
  dataclassobject@x <- data[sapply(data, is.numeric)]
  dataclassobject@y <- factor(data[sapply(data, is.factor)][,1])
  dataclassobject <- validatedata(dataclassobject)
  return(dataclassobject)
}




