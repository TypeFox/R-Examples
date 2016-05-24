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

#' Merge two objects of class \code{ProbDup}.
#' 
#' \code{MergeProbDup} merges two objects of class \code{ProbDup} into a single
#' one.
#' 
#' @param pdup1 An object of class \code{ProbDup}.
#' @param pdup2 An object of class \code{ProbDup}.
#' @return An object of class \code{ProbDup} with the merged list of  fuzzy, phonetic
#'   and semantic probable duplicate sets.
#' @seealso \code{\link[PGRdup]{ProbDup}}, \code{\link[PGRdup]{SplitProbDup}}
#' @examples
#' \dontrun{
#' #' # Load PGR passport database
#' GN <- GN1000
#' 
#' # Specify as a vector the database fields to be used
#' GNfields <- c("NationalID", "CollNo", "DonorID", "OtherID1", "OtherID2")
#' 
#' # Clean the data
#' GN[GNfields] <- lapply(GN[GNfields], function(x) DataClean(x))
#' y1 <- list(c("Gujarat", "Dwarf"), c("Castle", "Cary"), c("Small", "Japan"),
#' c("Big", "Japan"), c("Mani", "Blanco"), c("Uganda", "Erect"),
#' c("Mota", "Company"))
#' y2 <- c("Dark", "Light", "Small", "Improved", "Punjab", "SAM")
#' y3 <- c("Local", "Bold", "Cary", "Mutant", "Runner", "Giant", "No.",
#'         "Bunch", "Peanut")
#' GN[GNfields] <- lapply(GN[GNfields], function(x) MergeKW(x, y1, delim = c("space", "dash")))
#' GN[GNfields] <- lapply(GN[GNfields], function(x) MergePrefix(x, y2, delim = c("space", "dash")))
#' GN[GNfields] <- lapply(GN[GNfields], function(x) MergeSuffix(x, y3, delim = c("space", "dash")))
#' 
#' # Generate KWIC index
#' GNKWIC <- KWIC(GN, GNfields)
#' 
#' # Specify the exceptions as a vector
#' exep <- c("A", "B", "BIG", "BOLD", "BUNCH", "C", "COMPANY", "CULTURE", 
#'          "DARK", "E", "EARLY", "EC", "ERECT", "EXOTIC", "FLESH", "GROUNDNUT", 
#'          "GUTHUKAI", "IMPROVED", "K", "KUTHUKADAL", "KUTHUKAI", "LARGE", 
#'          "LIGHT", "LOCAL", "OF", "OVERO", "P", "PEANUT", "PURPLE", "R", 
#'          "RED", "RUNNER", "S1", "SAM", "SMALL", "SPANISH", "TAN", "TYPE", 
#'          "U", "VALENCIA", "VIRGINIA", "WHITE")
#'           
#' # Specify the synsets as a list
#' syn <- list(c("CHANDRA", "AH114"), c("TG1", "VIKRAM"))
#' 
#' # Fetch probable duplicate sets
#' GNdup <- ProbDup(kwic1 = GNKWIC, method = "a", excep = exep, fuzzy = TRUE,
#'                  phonetic = TRUE, encoding = "primary", 
#'                  semantic = TRUE, syn = syn)
#'                  
#' # Split the probable duplicate sets
#' GNdupSplit <- SplitProbDup(GNdup, splitat = c(10, 10, 10))
#' 
#' # Merge the split sets
#' GNdupMerged <- MergeProbDup(GNdupSplit[[1]], GNdupSplit[[3]])
#' 
#' }
#' @import data.table
#' @importFrom methods is
#' @export

MergeProbDup <- function(pdup1, pdup2) {

  if (!is(pdup1, "ProbDup")) {
    stop('"pdup1" is not of class ProbDup')
  }
  if (!is(pdup2, "ProbDup")) {
    stop('"pdup2" is not of class ProbDup')
  }
  if (!identical(attributes(pdup1), attributes(pdup2))) {
    stop('Attributes of "pdup1" and "pdup2" are not identical')
  }
  N <- length(seq_along(pdup1))
  outnames <- attributes(pdup1)$names
  out <- vector("list", length(outnames))
  names(out) <- outnames
  rm(outnames)
  attr(out, "method") <- attributes(pdup1)$method
  attr(out, "fields") <- attributes(pdup1)$fields
  class(out) <- class(pdup1)
  pdup1 <- lapply(pdup1, as.data.table)
  pdup2 <- lapply(pdup2, as.data.table)
  for (i in 1:N) {
    out[[i]] <- rbind.data.frame(pdup1[[i]], pdup2[[i]])
    out[[i]] <- unique(out[[i]], by = c("TYPE", "ID", "IDKW", "COUNT"))
    out[[i]] <- setDF(out[[i]])
    if (dim(out[[i]])[1] == 0) {
      out[i] <- list(NULL)
    }
  }
  return(out)
}
