### This file is part of 'PGRdup' package for R.

### Copyright (C) 2014, ICAR-NBPGR.
#
# PGRdup is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# PGRdup is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



#' Clean PGR passport data.
#' 
#' \code{DataClean} cleans the data in a character vector according to the 
#' conditions in the arguments.
#' 
#' This function aids in standardization and preparation of the PGR passport 
#' data for creation of a KWIC index with \code{\link[PGRdup]{KWIC}} function 
#' and the identification of probable duplicate accessions by the 
#' \code{\link[PGRdup]{ProbDup}} function. It cleans the character strings in 
#' passport data fields(columns) specified as the input character vector 
#' \code{x} according to the conditions in the arguments in the same order. If 
#' the input vector \code{x} is not of type character, it is coerced to a 
#' character vector.
#' 
#' This function is designed particularly for use with fields corresponding to 
#' accession names such as accession ids, collection numbers, accession names 
#' etc. It is essentially a wrapper around the \code{\link[base]{gsub}} base 
#' function with \code{\link[base]{regex}} arguments. It also converts all 
#' strings to upper case and removes leading and trailing spaces.
#' 
#' Commas, semicolons and colons which are sometimes used to separate multiple 
#' strings or names within the same field can be replaced with a single space 
#' using the logical arguments \code{fix.comma}, \code{fix.semcol} and 
#' \code{fix.col} respectively.
#' 
#' Similarly the logical argument \code{fix.bracket} can be used to replace all 
#' brackets including parenthesis, square brackets and curly brackets with
#' space.
#' 
#' The logical argument \code{fix.punct} can be used to remove all punctuation 
#' from the data.
#' 
#' \code{fix.space} can be used to convert all space characters such as tab, 
#' newline, vertical tab, form feed and carriage return to spaces and finally 
#' convert multiple spaces to single space.
#' 
#' \code{fix.sep} can be used to merge together accession identifiers 
#' composed of alphabetic characters separated from as series of digits by a 
#' space character. For example IR 64, PUSA 256 etc.
#' 
#' \code{fix.leadzero} can be used to remove leading zeros from accession name 
#' fields to facilitate matching to identify probable duplicates. e.g. IR0064 -> 
#' IR64
#' 
#' @param x A character vector. If not, coerced to character by 
#'   \code{as.character}.
#' @param fix.comma logical. If \code{TRUE}, all the commas are replaced by 
#'   space (see \strong{Details}).
#' @param fix.semcol logical. If \code{TRUE}, all the semicolons are replaced by
#'   space (see \strong{Details}).
#' @param fix.col logical. If \code{TRUE}, all the colons are replaced by space 
#'   (see \strong{Details}).
#' @param fix.bracket logical. If \code{TRUE}, all the brackets are replaced by 
#'   space (see \strong{Details}).
#' @param fix.punct logical. If \code{TRUE}, all punctuation characters are 
#'   removed (see \strong{Details}).
#' @param fix.space logical. If \code{TRUE}, all space characters are replaced 
#'   by space and multiple spaces are converted to single space (see 
#'   \strong{Details}).
#' @param fix.sep logical. If \code{TRUE}, space between alphabetic characters 
#'   followed by digits is removed (see \strong{Details}).
#' @param fix.leadzero logical. If \code{TRUE}, leading zeros are removed (see 
#'   \strong{Details}).
#' @return A character vector with the cleaned data converted to upper case. 
#'   \code{NAs} if any are converted to blank strings.
#' @seealso \code{\link[base]{gsub}}, \code{\link[base]{regex}}, 
#'   \code{\link[PGRdup]{MergeKW}}, \code{\link[PGRdup]{KWIC}}, 
#'   \code{\link[PGRdup]{ProbDup}}
#' @examples
#' names <- c("S7-12-6", "ICG-3505", "U 4-47-18;EC 21127", "AH 6481", "RS   1",
#'            "AK 12-24", "2-5 (NRCG-4053)", "T78, Mwitunde", "ICG 3410",
#'            "#648-4 (Gwalior)", "TG4;U/4/47/13", "EC0021003")
#' DataClean(names)
#' @export
DataClean <- function(x, fix.comma = TRUE, fix.semcol = TRUE, fix.col = TRUE,
                       fix.bracket = TRUE, fix.punct = TRUE, fix.space = TRUE,
                       fix.sep = TRUE, fix.leadzero = TRUE) {
  if (!is.character(x)) {
    warning("x is not of type character; coerced to character")
    x <- as.character(x)
  }
    x[is.na(x)]   <- ""  # Convert NAs to empty strings
    # Convert all strings to upper case
    x <- vapply(x, FUN = toupper, FUN.VALUE = "character")
    if (fix.comma) {
    x <- gsub(pattern = ",", replacement = " ", x)  # Replace "," by space
    }
    if (fix.semcol) {
    x <- gsub(pattern = ";", replacement = " ", x)  # Replace ";" by space
    }
    if (fix.col) {
    x <- gsub(pattern = ":", replacement = " ", x)  # Replace ":" by space
    }
    # Replace brackets by space
    if (fix.bracket) {
    x <- gsub(pattern = "[](){}[]", replacement = " ", x)
    }
    if (fix.punct) {
    # Remove all punctuations
    x <- gsub(pattern = "[[:punct:]]", replacement = "", x, perl = TRUE)
    }
    if (fix.space) {
    # replace all space characters with space (" ")
    x <- gsub(pattern = "[[:space:]]", replacement = " ", x)
    # replace multiple spaces with single space
    x <- gsub(pattern = "([[:space:]])\\1+", replacement = "\\1", x)
    }
    if (fix.sep) {
    # Remove space when alphabetic characters followed by digits
    x <- gsub(pattern = "([a-zA-Z]) ([0-9])", replacement = "\\1\\2", x,
              perl = TRUE)
    }
    if (fix.leadzero) {
    # Remove leading zeros
    x <- gsub("(?<![0-9])0+", "", x, perl = TRUE)
    }
    # Remove leading and trailing spaces if any.
    x <- gsub("^\\s+|\\s+$", "", x)
    # Remove names attribute
    names(x) <- NULL
return(x)
}
