#' Fill missing values in data set (or matrix) columns for which there is a
#' single unique value
#' 

#' 
#' Looks for columns in a data set that have a single unique non-missing value
#' and fills in all \code{NA} with that value
#' 
#' 
#' @param data data.frame
#' @return data.frame
#' @author Devin S. Johnson
#' @examples
#' 
#' #library(crawl)
#' data1 <- data.frame(constVals=rep(c(1,NA),5), vals=1:10)
#' data1[5,2] <- NA
#' data1
#' data2 <- fillCols(data1)
#' data2
#' 
#' mat1 <- matrix(c(rep(c(1,NA),5), 1:10), ncol=2)
#' mat1[5,2] <- NA
#' mat1
#' mat2 <- fillCols(mat1)
#' mat2
#' @export
"fillCols" <- function(data) {
   nc <- ncol(data)
   getConst <- function(vec) {
      vals <- unique(vec)
      return(length(vals[!is.na(vals)])==1)
   }
   constCol <- apply(data, 2, getConst)
   data[,constCol] <- data[1,constCol]
   return(data)
}
