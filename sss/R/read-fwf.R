# 
# Author: Andrie
#----------------------------------------------------------------------------------



#' Read fixed-width files quickly.
#' 
#' Experimental replacement for read.fwf that runs much faster.  However, it is much less flexible than read.fwf.
#' 
#' @param file Character vector: name of file 
#' @param widths Numeric vector: column widths. Negative numbers mean "skip this many columns". Use an NA as the final element if there are likely to be extra characters at the end of each row after the last one that you're interested in. 
#' @param col.names names for the columns that are NOT skipped 
#' @param colClasses can be used to control type conversion; see \code{\link{read.table}}. It is an optional vector whose names must be part of col.names. There is one extension of the \code{\link{read.table}} rules:a colClasses string starting \code{\link{POSIXct}}. will trigger automatic conversion to POSIXct, using the rest of the string as the format specifier. 
#' @param tz used in auto-conversion to \code{\link{POSIXct}} when \code{colClasses} is set 
#' @param dec the character to be assumed for decimal points. Passed to \code{\link[utils]{type.convert}}
#' @param ... ignored 
fast.read.fwf <- function (file, widths, col.names = NULL, colClasses = NA,  
                           tz = "", dec = ".", ...){
  
  if(is.null(col.names)) col.names <- paste("V", seq_along(widths), sep = "")
  
  fs <- file.size(file)
  
  ll <- readLines(file, encoding = "UTF-8")
  nl <- length(ll)
  rw <- lapply(ll, charToRaw)
  
  fields <- vapply(rw, readChar, abs(widths), useBytes = TRUE, 
                   FUN.VALUE = character(length(widths)))
  # length(fields)
  
  fields <- matrix(fields, nrow = nl, ncol = length(widths), byrow = TRUE)
  #colnames(fields) <- col.names
  
  modColClass <- function(i){
    if (is.na(colClasses[i])) 
    type.convert(fields[, i], as.is = TRUE, dec = dec, na.strings = character(0L))
    else switch(colClasses[i],
                "factor" = as.factor(fields[, i]),
                "Date" = as.Date(fields[, i], format = "%Y%m%d"),
                "POSIXct" = as.POSIXct(fields[, i]),
                "logical" = as.logical(as.numeric(fields[, i])),
                suppressWarnings(
                  methods::as(fields[, i], colClasses[i])
                )
                
    )
  }
  
  df2 <- lapply(seq_len(ncol(fields)), modColClass)
  df2 <- fastdf(df2)
  names(df2) <- col.names
  df2
}
