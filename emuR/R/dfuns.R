##' Find the time and position of a data element.
##' 
##' Finds the time and position of a data element.
##' 
##' The dataset returned from \code{track} or \code{frames} consists of a
##' matrix of data (the \code{data} component) and two index components
##' (\code{index} and \code{ftime}). The data for all segments is concatenated
##' together in \code{$data}.  This function can be used to find out which
##' segment a particular row of \code{$data} corresponds to.
##' 
##' @param dataset A dataset returned by \code{track} or \code{frames}.
##' @param datanum An integer, an index into the \code{data} component of
##' \code{dataset}.
##' @return The segment number which contains the element \code{datanum} of
##' \code{dataset$data}.
##' @seealso track, frames
##' @keywords misc
##' @export frames.time
"frames.time" <- function(dataset, datanum)
{
  ## return the time and the number of the segment element
  ## that the datanum refers to
  if(is.matrix(dataset$ftime) == FALSE) dataset$ftime <- rbind(dataset$ftime)
  if(is.matrix(dataset$index) == FALSE)  dataset$index <- rbind(dataset$index)
  nums <- seq(1, nrow(dataset$ftime))
  incl <- dataset$index[, 1] <= datanum & dataset$index[, 2] >= datanum
  retv <- NULL
  segnum <- nums[incl]
  percent <- (datanum - dataset$index[segnum, 1])/
             (dataset$index[segnum, 2] - dataset$index[segnum, 1])
  retv$segnum <- segnum
  retv$time <- dataset$ftime[segnum, 1] + 
               percent * (dataset$ftime[segnum, 2] - dataset$ftime[segnum, 1])
  retv
}










##' Get data for a given time
##' 
##' Gets data for a given time
##' 
##' 
##' @param timeval A time in milliseconds
##' @param dataset A trackdata object as returned by \code{track}.
##' @return The element number of \code{trackdata$data} corresponding to
##' \code{time}
##' @seealso track, frames
##' @keywords misc
##' @export get.time.element
"get.time.element"<- function(timeval, dataset)
{
  ## timeval: a time in milliseconds
  ## dataset: a data structure consisting of $data, $ftime, $index
  ## returns the element number of dataset$data corresponding to timeval
  numrows <- nrow(dataset$ftime)
  left <- dataset$ftime[1, 1]
  right <- dataset$ftime[numrows, 2]
  left.i <- dataset$index[1, 1]
  right.i <- dataset$index[numrows, 2]
  round(((timeval - left)/(right - left)) * (right.i - left.i)) + 1
}
