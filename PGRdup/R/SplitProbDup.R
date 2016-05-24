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

#' Split an object of class \code{ProbDup}.
#' 
#' \code{SplitProbDup} splits an object of class \code{ProbDup} into two on the 
#' basis of set counts.
#' 
#' @param pdup An object of class \code{ProbDup}.
#' @param splitat A vector of 3 integers indicating the set count at which 
#'   Fuzzy, Phonetic and Semantic duplicate sets in \code{pdup} are to be split.
#' @return A list with the the divided objects of class \code{ProbDup} 
#'   (\code{pdup1} and \code{pdup2}) along with the corrsponding lists of
#'   accessions present in each (\code{list1} and \code{list2}).
#' @seealso \code{\link[PGRdup]{ProbDup}}, \code{\link[PGRdup]{MergeProbDup}}
#' @examples
#' \dontrun{
#' # Load PGR passport database
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
#' }
#' @import data.table
#' @importFrom methods is
#' @export


SplitProbDup <- function(pdup, splitat = c(30, 30, 30)) {
  if (!is(pdup, "ProbDup")) {
    stop('"pdup" is not of class ProbDup')
  }
  N <- length(seq_along(pdup))
  outnames <- attributes(pdup)$names
  out1 <- vector("list", length(outnames))
  names(out1) <- outnames
  rm(outnames)
  attr(out1, "method") <- attributes(pdup)$method
  attr(out1, "fields") <- attributes(pdup)$fields
  class(out1) <- class(pdup)
  out2 <- out1
  pdup <- lapply(pdup, function(x) if (!is.null(x)) as.data.table(x))
  list1 <- data.table(TYPE = vector(mode = "character"),
                      V1 = vector(mode = "character"))
  list2 <- data.table(TYPE = vector(mode = "character"),
                      V1 = vector(mode = "character"))
  for (i in 1:N) {
    if (!is.null(pdup[[i]])) {
      out1[[i]] <- subset(pdup[[i]], COUNT <= splitat[[i]])
      out2[[i]] <- subset(pdup[[i]], COUNT > splitat[[i]])
      out1[[i]][, SET_NO := as.numeric(as.factor(SET_NO))]
      out2[[i]][, SET_NO := as.numeric(as.factor(SET_NO))]
      if (!dim(out1[[i]])[1] == 0) {
        list1 <- rbind(list1,
                       out1[[i]][, .(unlist(strsplit(ID, ", ", TRUE))),
                                 by = list(TYPE)])
      }
      if (!dim(out2[[i]])[1] == 0) {
        list2 <- rbind(list2,
                       out2[[i]][, .(unlist(strsplit(ID, ", ", TRUE))),
                                 by = list(TYPE)])
      }
      out1[[i]] <- setDF(out1[[i]])
      out2[[i]] <- setDF(out2[[i]])
      if (dim(out1[[i]])[1] == 0) {
        out1[i] <- list(NULL)
        }
      if (dim(out2[[i]])[1] == 0) {
        out2[i] <- list(NULL)
        }
    }
  }
  if (dim(list1)[1] == 0) {
    list1 <- NULL
    }
  if (dim(list2)[1] == 0) {
    list2 <- NULL
    }
  if (!is.null(list1)) {
    list1 <- unique(list1, by = colnames(list1))
    list1[, PRIM_ID := gsub("\\[K1\\]", "", V1, perl = TRUE)]
    list1[, PRIM_ID := gsub("\\[K2\\]", "", PRIM_ID, perl = TRUE)]
    list1[, V1 := gsub("(\\]).*", "\\1", V1, perl = TRUE)]
    setnames(list1, "V1", "K")
    setDF(list1)
  }
  if (!is.null(list2)) {
    list2 <- unique(list2, by = colnames(list2))
    list2[, PRIM_ID := gsub("\\[K1\\]", "", V1, perl = TRUE)]
    list2[, PRIM_ID := gsub("\\[K2\\]", "", PRIM_ID, perl = TRUE)]
    list2[, V1 := gsub("(\\]).*", "\\1", V1, perl = TRUE)]
    setnames(list2, "V1", "K")
    setDF(list2)
  }
  out <- list(pdup1 = NULL, list1 = NULL,
              pdup2 = NULL, list2 = NULL)
  if (is.element(FALSE, lapply(out1, is.null))) {
    out[[1]] <- out1
    }
  if (!is.null(list1)) {
    out[[2]] <- list1
    }
  if (is.element(FALSE, lapply(out2, is.null))) {
    out[[3]] <- out2
    }
  if (!is.null(list2)) {
    out[[4]] <- list2
    }
  rm(out1, out2, list1, list2)
  return(out)
}
