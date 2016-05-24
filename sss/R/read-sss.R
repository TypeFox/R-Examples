#' Reads a triple-s XML (asc) data file, as specified by the triple-s XML standard.
#'
#' This function reads and parses a .sss XML metadata file as well as its associated .asc data file. The .sss standard defines a standard survey structure
#'
#' @param sssFilename Character string: name of .sss file containing the survey metadata
#' @param ascFilename Character string: name of .asc file containing survey data
#' @param sep Character vector defining the string that separates question and subquestion labels, e.g. \code{c("Q_1", "Q_2")}
#' @return
#' A data frame with one element (column) for each variable in the data set.
#' The data.frame contains several attributes:
#' 
#' \describe{
#' \item{variable.labels}{a named list of value labels with one element per variable, either NULL or a names character vector}
#' }
#' @keywords read
#' @references http://www.triple-s.org/
#' @export 
#' @examples
#' # Not executed
#' # read.sss("sample.sss, sample.asc")
read.sss <- function(sssFilename, ascFilename, sep = "_"){
  message("Reading SSS metadata")
  switch(class(sssFilename),
         "character" = {
           if(!file.exists(sssFilename)) stop("File doesn't exist: ", sssFilename)
           if(!file.exists(ascFilename)) stop("File doesn't exist: ", ascFilename)
           doc <- readSSSmetadata(sssFilename)
           sss <- parseSSSmetadata(doc)
         }, 
         "XMLDocumentContent" = {
           sss <- parseSSSmetadata(sssFilename)
         }, 
         stop("SSSfilename not recognised as either a file or an XML object")
  )
  
  sss$variables <- splitSSS(sss$variable, sep)
  
  message("Reading SSS data")
  #asc <- readSSSdata(ascFilename)
  
  ascWidth <- sss$variables$colWidth
  
  types <- c(single = "character",
              multiple = "character",
              character = "character", 
              logical = "logical",
              numeric = "numeric", 
              quantity = "numeric",
              date = "Date"
              )
  ascType <- types[sss$variables$type]
  ascType[sss$variables$type == "multiple"] <- "numeric"
  ascType[sss$variables$type == "multiple" & sss$variables$subfields > 0] <- "character"
  
  ascNames <- sss$variables$name
  
  dat <- switch(sss$format, 
                csv = 
                  read.csv(file = ascFilename,
                           skip = sss$skip,
                           header = FALSE,
                           col.names = ascNames,
                           colClasses = "character",
                           stringsAsFactors = FALSE
                  ),
                fixed = 
                  fast.read.fwf(file = ascFilename, 
                                widths = ascWidth, 
                                colClasses = ascType, 
                                col.names = ascNames
                  )
                
  )
  dat <- changeValues(sss, dat)
  dat <- addQtext(sss, dat)
  dat
}

