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



#' Merge keyword strings
#' 
#' These functions merge keyword strings separated by delimiters such as space, 
#' period or dash in a character vector into single keyword strings.
#' 
#' These functions aid in standardization of relevant data fields(columns) in 
#' PGR passport data for creation of a KWIC index with 
#' \code{\link[PGRdup]{KWIC}} function and subsequent identification of probable
#' duplicate accessions by the \code{\link[PGRdup]{ProbDup}} function.
#' 
#' It is recommended to run this function before using the 
#' \code{\link[PGRdup]{DataClean}} function on the relevant data fields(columns)
#' of PGR passport databases.
#' 
#' \code{MergeKW} merges together pairs of strings specified as a list in 
#' argument \code{y} wherever they exist in a character vector. The second 
#' string in the pair is merged even when it is followed by a number.
#' 
#' \code{MergePrefix} merges prefix strings specified as a character vector in 
#' argument \code{y} to the succeeding root word, wherever they exist in a 
#' character vector.
#' 
#' \code{MergeSuffix} merges suffix strings specified as a character vector in 
#' argument \code{y} to the preceding root word, wherever they exist in a 
#' character vector. The suffix strings which are followed by numbers are also 
#' merged.
#' 
#' @param x A character vector. If not, coerced to character by 
#'   \code{as.character}.
#' @param y A list of character vectors with pairs of strings that are to be 
#'   merged (for \code{MergeKW}) or a character vector of strings which are to 
#'   be merged to succeeding string (for \code{MergePrefix}) or the preceding 
#'   string (for \code{MergeSuffix}). If not of type character, coerced by 
#'   \code{as.character}.
#' @param delim Delimiting characters to be removed between keywords.
#' @return A character vector of the same length as \code{x} with the required 
#'   keyword strings merged.
#' @examples
#' names <- c("Punjab Bold", "Gujarat- Dwarf", "Nagpur.local", "SAM COL 144",
#'            "SAM COL--280", "NIZAMABAD-LOCAL", "Dark Green Mutant",
#'            "Dixie-Giant", "Georgia- Bunch", "Uganda-erect", "Small Japan",
#'            "Castle  Cary", "Punjab erect", "Improved small japan",
#'            "Dark Purple")
#'
#' # Merge pairs of strings
#' y1 <- list(c("Gujarat", "Dwarf"), c("Castle", "Cary"), c("Small", "Japan"),
#'            c("Big", "Japan"), c("Mani", "Blanco"), c("Uganda", "Erect"),
#'            c("Mota", "Company"))
#' names <- MergeKW(names, y1, delim = c("space", "dash", "period"))
#'
#' # Merge prefix strings
#' y2 <- c("Light", "Small", "Improved", "Punjab", "SAM")
#' names <- MergePrefix(names, y2, delim = c("space", "dash", "period"))
#'
#' # Merge suffix strings
#' y3 <- c("Local", "Bold", "Cary", "Mutant", "Runner", "Giant", "No.",
#'         "Bunch", "Peanut")
#' names <- MergeSuffix(names, y3, delim = c("space", "dash", "period"))
#' @seealso \code{\link[PGRdup]{DataClean}}, \code{\link[PGRdup]{KWIC}}, 
#'   \code{\link[PGRdup]{ProbDup}}
#' @name MergeKW

#' @rdname MergeKW
#' @export
MergeKW <- function(x, y, delim = c("space", "dash", "period")) {
  if (!is.character(x)) {
    warning("x is not of type character; coerced to character")
    x <- as.character(x)
  }
  if (!is.list(y)) {
    stop("y is not a list")
  }
  if (is.element(FALSE, as.logical(lapply(y, function(x) length(x) == 2)))) {
    stop("list y is not in the appropriate format")
  }
  if (is.element(FALSE, as.logical(lapply(y, function(x) is.character(x))))) {
    warning("list y had non character vectors; coerced to character")
    y <-  as.logical(lapply(y, function(x) as.character(x)))
  }
  delim <- match.arg(delim, c("space", "dash", "period"), several.ok = TRUE)
  # replace all space characters with space (" ")
  x <- gsub(pattern = "[[:space:]]", replacement = " ", x)
  y <- unique(y)
  # Escape all Regex special characters in y
  y <- lapply(y, function(x) gsub(pattern = "([.|()\\^{}+$*?]|\\[|\\])",
                                  replacement = "\\\\\\1", x))
  options <- c("\\s", "-", ".")
  options2 <- logical(length = 3)
  if (is.element("space", delim)) {
    options2[1] <- TRUE
  }
  if (is.element("dash", delim)) {
    options2[2] <- TRUE
  }
  if (is.element("period", delim)) {
    options2[3] <- TRUE
  }
  p <-  paste0(options[options2], collapse = "")
  sapply(seq_len(length(y)), function(i){
    pat <- paste0("(?i)(?<=^", y[[i]][1], ")[", p, "]+(?=", y[[i]][2],
                  ")[[:digit:]]*\\b|\\b(?<=", y[[i]][1], ")[", p,
                  "]+(?=",  y[[i]][2], "[[:digit:]]*$)", collapse = "")
    x <<- gsub(pat, "", x, perl = TRUE)
    }
  )
  return(x)
}

#' @rdname MergeKW
#' @export
MergePrefix <- function(x, y, delim = c("space", "dash", "period")) {
  if (!is.character(x)) {
    warning("x is not of type character; coerced to character")
    x <- as.character(x)
  }
  if (!is.character(y)) {
    warning("x is not of type character; coerced to character")
    y <- as.character(y)
  }
  delim <- match.arg(delim, c("space", "dash", "period"), several.ok = TRUE)
  # replace all space characters with space (" ")
  x <- gsub(pattern = "[[:space:]]", replacement = " ", x)
  # Escape all Regex special characters in y
  y <- unique(toupper(y))
  y <- gsub(pattern = "([.|()\\^{}+$*?]|\\[|\\])",
            replacement = "\\\\\\1", y)
  options <- c("\\s", "-", ".")
  options2 <- logical(length = 3)
  if (is.element("space", delim)) {
    options2[1] <- TRUE
  }
  if (is.element("dash", delim)) {
    options2[2] <- TRUE
  }
  if (is.element("period", delim)) {
    options2[3] <- TRUE
  }
  p <-  paste0(options[options2], collapse = "")
  sapply(seq_len(length(y)), function(i){
    x <<- gsub(paste0("(?i)\\b(?<=^", y[i],
                      ")[", p, "]+"), "", x, perl = TRUE)
    }
  )
  return(x)
}

#' @rdname MergeKW
#' @export
MergeSuffix <- function(x, y, delim = c("space", "dash", "period")) {
  if (!is.character(x)) {
    warning("x is not of type character; coerced to character")
    x <- as.character(x)
  }
  if (!is.character(y)) {
    warning("x is not of type character; coerced to character")
    y <- as.character(y)
  }
  delim <- match.arg(delim, c("space", "dash", "period"), several.ok = TRUE)
  # replace all space characters with space (" ")
  x <- gsub(pattern = "[[:space:]]", replacement = " ", x)
  # Escape all Regex special characters in y
  y <- unique(toupper(y))
  y <- gsub(pattern = "([.|()\\^{}+$*?]|\\[|\\])",
            replacement = "\\\\\\1", y)
  options <- c("\\s", "-", ".")
  options2 <- logical(length = 3)
  if (is.element("space", delim)) {
    options2[1] <- TRUE
  }
  if (is.element("dash", delim)) {
    options2[2] <- TRUE
  }
  if (is.element("period", delim)) {
    options2[3] <- TRUE
  }
  p <-  paste0(options[options2], collapse = "")
  sapply(seq_len(length(y)), function(i){
      x <<- gsub(paste0("(?i)[", p, "]+(?=", y[i], "[[:digit:]]*",
                        "\\b)"), "", x, perl = TRUE)
      }
    )
  return(x)
}
