# Copyright (C) 2010-2014 Leo Lahti, Juuso Parkkinen, Joona Lehtomaki 
# <ropengov.github.com>. All rights reserved.

# This program is open source software; you can redistribute it and/or modify 
# it under the terms of the FreeBSD License (keep this notice): 
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' List the open data files available from Statistics Finland
#' 
#' Arguments:
#'  @param format "px", "csv", "xml" depending on the desired format of the data files
#'
#' Returns:
#'  @return table
#'
#' @export
#' @references
#' See citation("statfi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples df <- list_statfi_files()
#' @keywords utilities

list_statfi_files <- function (format = "px") {

  if (format == "px") {
    url <- "http://pxweb2.stat.fi/database/StatFin/StatFin_rap.csv"
  } else {  
    url <- paste("http://pxweb2.stat.fi/database/StatFin/StatFin_rap_", format, ".csv", sep = "")
  }

  message(paste("Downloading", url))
  d <- try(read.csv(url, sep = ";", encoding = "latin1"))
  if (length(grep("^Error", d)) == 1) {warning(paste("No files available at", url)); return(character(0))}

  tab <- apply(d, 2, function (x) {iconv(x, from = "latin1", to = "utf-8")})
  tab <- as.data.frame(tab)

  tab$File <- as.character(tab$File)
  tab$size <- as.numeric(as.character(tab$size)) 
  tab$created <- as.character(tab$created) 
  tab$updated <- as.character(tab$updated) 
  tab$tablesize <- as.character(tab$tablesize) 
  tab$TITLE <- as.character(tab$TITLE)
  tab$DESCRIPTION <- as.character(tab$DESCRIPTION)

  tab

}


#' List the open data files available from Eurostat
#' 
#' Arguments:
#'  @param ... Arguments to be passed
#'
#' Returns:
#'  @return table
#'
#' @export
#' @references
#' See citation("statfi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples # df <- list_eurostat_files()
#' @keywords utilities

list_eurostat_files <- function (...) {

  url <- "http://pxweb2.stat.fi/database/StatFin/StatFin_rap.csv"

  message(paste("Downloading", url))
  d <- try(read.csv(url, sep = ";", encoding = "latin1"))
  if (length(grep("^Error", d)) == 1) {warning(paste("No files available at", url)); return(character(0))}

  tab <- apply(d, 2, function (x) {iconv(x, from = "latin1", to = "utf-8")})
  tab <- as.data.frame(tab)

  tab$File <- as.character(tab$File)
  tab$size <- as.numeric(as.character(tab$size)) 
  tab$created <- as.character(tab$created) 
  tab$updated <- as.character(tab$updated) 
  tab$variables <- as.numeric(tab$variables)
  tab$tablesize <- as.character(tab$tablesize) 
  tab$TITLE <- as.character(tab$TITLE)
  tab$DESCRIPTION <- as.character(tab$DESCRIPTION)

  tab

}




#' Get PC Axis data with custom preprocessing for PC Axis 
#' files from Statistics Finland (Tilastokeskus) http://www.stat.fi/
#'
#' Arguments:
#'  @param url or local file name of the StatFi file
#'  @param format One of the following: "px", "csv", "xml". Specifies the desired format of the source file.
#'  @param verbose verbose
#'
#' Returns:
#'  @return data.frame
#'
#' @details If the reading of PX file fails, CSV is used instead.
#'
#' @export
#' @references
#' See citation("statfi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples \dontrun{px <- get_statfi("http://pxweb2.stat.fi/database/StatFin/vrm/synt/080_synt_tau_203.px")}
#' @keywords utilities

get_statfi <- function (url, format = "px", verbose = TRUE) {

  if (format == "px") {

    url <- gsub("\\.csv", "\\.px", url)
    url <- gsub("\\.xml", "\\.px", url)

    # If URL is given, read the data into PX object
    if (is_url(url)) {
      message(paste("Reading StatFi data from ", url))
      px <- read_px(url)
    }

    # Convert to data.frame 
    if (class(px) == "px") { 
      df <- as.data.frame(px) 
    }

  } else if (format == "xml") {

    warning("xml not yet implemented for statfi; using csv instead")

    url <- gsub("\\.px", "\\.csv", url)
    url <- gsub("\\.xml", "\\.csv", url)
    df <- read.csv(url, encoding = "latin1", as.is = T, colClasses = 'character', sep = ";"); 

  } else if (format == "csv") {

    # TODO
    url <- gsub("\\.px", "\\.csv", url)
    url <- gsub("\\.xml", "\\.csv", url)
    df <- read.csv(url, encoding = "latin1", as.is = T, colClasses = 'character', sep = ";"); 

  }

  # Some preprocessing for field names
  # TODO Improve! Also convert all to UTF8
  fields <- c("Alue", "Kunta")
  for (nam in intersect(fields, colnames(df))) {
    df[[nam]] <- sapply(df[[nam]], function (x) {strsplit(as.character(x), " - ")[[1]][[1]]})
  }    
  fields <- c("Vuosi")
  for (nam in intersect(fields, colnames(df))) {
    df[[nam]] <- as.numeric(as.character(df[[nam]]))
  }    

  df

}



