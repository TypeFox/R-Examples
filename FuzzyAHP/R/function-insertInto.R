#' Function that inserts FuzzyData on specific column number in another FuzzyData
#'
#' @description
#' This methods construct object \code{\linkS4class{FuzzyData}} based on two \code{\linkS4class{FuzzyData}} and index.
#' The functions merges the sources into single output. It inserts data2 onto specific position (defined by index) in
#' data1.
#'
#' @param data1 An object of \code{\linkS4class{FuzzyData}} - data to insert into.
#' @param data2 An object of \code{\linkS4class{FuzzyData}} - data that should be inserted.
#' @param index An integer specifying the column position on which data2 should be inserted.
#'
#' @return An object of class \code{\linkS4class{FuzzyData}}
#'
#' @export
#' @rdname insertInto-methods
#' @name insertInto
setGeneric("insertInto",
           function(data1, data2, index) standardGeneric("insertInto"))

#' @rdname insertInto-methods
#' @aliases insertInto,FuzzyData,FuzzyData,integer-method
#'
#' @examples
#' values = (c(1,2,3,4,5,6,7,8,9))
#' values = matrix(values, nrow = 3, ncol = 3, byrow = TRUE)
#' fData = fuzzyData(values, single.value = FALSE)
#'
#' v = (c(15,16,17))
#' v = matrix(v, nrow = 3, ncol = 1, byrow = TRUE)
#' fd = fuzzyData(v, single.value = FALSE)
#'
#' fData = insertInto(fData, fd, 2)
#'
setMethod(
  f= "insertInto",
  signature(data1 = "FuzzyData", data2 = "FuzzyData"),
  definition=function(data1, data2, index)
  {
    if(index != as.integer(index)){
      stop("Index needs to be an integer value.")
    }

    if(index > ncol(data1@fnModal)){
      stop(paste("Index outside of range. Can not put data2 on position ", index, " when there are only ", ncol(data1@fnModal),
                 " of colums in data1.", sep = ""))
    }

    if(ncol(data2@fnModal)>1){
      stop("Can only insert data2 into specific position if it contains just one fuzzy value.")
    }

    if(ncol(data1@fnModal)==1){
      if(index == 1){
        return(bindColums(data2, data1))
      }
    }

    sequence = seq(ncol(data1@fnModal), index, -1)

    data1@fnMin = cbind(data1@fnMin, 0)
    data1@fnModal = cbind(data1@fnModal, 0)
    data1@fnMax = cbind(data1@fnMax, 0)

    for(i in sequence){
      data1@fnMin[,i+1] = data1@fnMin[,i]
      data1@fnModal[,i+1] = data1@fnModal[,i]
      data1@fnMax[,i+1] = data1@fnMax[,i]
    }

    data1@fnMin[,index] = data2@fnMin[,1]
    data1@fnModal[,index] = data2@fnModal[,1]
    data1@fnMax[,index] = data2@fnMax[,1]

    return(data1)
  }
)
