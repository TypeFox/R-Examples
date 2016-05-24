#' Function to combine water-quality sample data and continuous (daily) 
#' ancillary variables and drop unnecessary columns.
#'
#' @name combineData
#' @title Combine water-quality sample data and continuous ancillary 
#' variables
#' @param qwdat is the dataset containing water-quality sample data with 
#' columns that begin with a P (or other user-defined indicator) followed 
#' by alphanumeric characters.  These columns are concentration data.  In 
#' addition there need to be columns that begin with an R (or other user-
#' defined indicator) followed by alphanumeric characters that match those 
#' of the associated concentration data.  The R columns contain data 
#' qualification codes.  See example datasets for more information about 
#' the data format, \link{IllRivValleyCty} and \link{qwMoRivOmaha}.
#' @param qwcols is a character vector with column headings for 
#' a station (location) identifier, a dates column identifier, beginning 
#' of column headers for remarks code (default is R), and beginning of 
#' column headers for concentration data (default is P for parameter).
#' @param cqwdat is the dataset containing variables that can be used as 
#' explanatory variables for the seawaveQ model.  See example dataset
#' for more information about the data format \link{cqwMoRivOmaha}.  These 
#' are daily values with no missing values allowed between the first and 
#' the last date in the dataset.
#' @return a dataframe
#' @format a dataframe with the number of rows equal to the number of rows
#' in the dataframe indicated by qwdat.  The number of columns depend 
#' on the two input data frames.  Minimally there will be a station 
#' identification column, a dates column, a column of qualification codes,
#' and a column of water-quality data.
#' @note The columns indicated by qwcols[1:2] are used to
#' combined the datasets.  The first column is the station identifier and 
#' the second column is the date column.  These two column headings must 
#' be the same in the two datasets being combined and the dates in the 
#' datasets being combined must be for the class Date and be in the same 
#' format.
#' @keywords manip
#' @author Karen R. Ryberg and Aldo V. Vecchia
#' @export
#' @examples
#' data(swData)
#' MoRivOmaha<-combineData(qwdat=qwMoRivOmaha, cqwdat=cqwMoRivOmaha, 
#' qwcols=c("staid", "dates", "R", "P"))

combineData <- function(qwdat, cqwdat, 
                        qwcols=c("staid", "dates", "R", "P")) {
  # combine concentration data with continous data and anomalies
  pattern <- paste("(", qwcols[1], ")|(", qwcols[2], ")|(", qwcols[3],"|", 
                   qwcols[4], ")[[:alnum:]]+", sep="")
  mycols <- grep(pattern, dimnames(qwdat)[[2]])
  # drop extract columns
  qwdat <- qwdat[,mycols]
  dat <- merge(qwdat, cqwdat, by.x=qwcols[1:2], by.y=qwcols[1:2], 
               all.x=TRUE )
  row.names(dat) <- NULL
  dat
}
