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



#' Reconstruct an object of class ProbDup
#' 
#' \code{ReconstructProbDup} reconstructs a data frame of probable duplicate 
#' sets created using the function \code{ReviewProbDup} and subjected to manual 
#' clerical review, back into an object of class \code{ProbDup}.
#' 
#' A data frame created using the function \code{\link[PGRdup]{ReviewProbDup}} 
#' from an object of class \code{ProbDup} for manual clerical review of 
#' identified probable duplicate sets can be reconstituted back to the same 
#' object after the review using this function. The instructions for modifying 
#' the sets entered in the appropriate format in the columns \code{DEL} and 
#' \code{SPLIT} during clerical review are taken into account for reconstituting
#' the probable duplicate sets.
#' 
#' Any records with \code{Y} in column \code{DEL} are deleted and records with 
#' identical integers in the column \code{SPLIT} other than the default \code{0}
#' are reassembled into a new set.
#' 
#' @param rev A data frame with the the core columns(fields) \code{SET_NO}, 
#'   \code{TYPE}, \code{K}, \code{PRIM_ID}, \code{DEL}, \code{SPLIT}, 
#'   \code{COUNT} and \code{IDKW}
#' @return An object of class \code{ProbDup} with the modified fuzzy,
#'   phonetic and semantic probable duplicate sets according to the instructions
#'   specified under clerical review.
#' @examples
#' \dontrun{
#'
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
#' # Get disjoint probable duplicate sets of each kind
#' disGNdup <- DisProbDup(GNdup, combine = NULL)
#' 
#' # Get the data frame for reviewing the duplicate sets identified
#' RevGNdup <- ReviewProbDup(pdup = disGNdup, db1 = GN1000,
#'                           extra.db1 = c("SourceCountry", "TransferYear"), 
#'                           max.count = 30, insert.blanks = TRUE)
#' # Examine and review the duplicate sets using edit function
#' RevGNdup <- edit(RevGNdup)
#' 
#' # Examine and make changes to a set
#' subset(RevGNdup, SET_NO==12 & TYPE=="P", select= c(IDKW, DEL, SPLIT))
#' RevGNdup[c(110, 112, 114, 118, 121, 122, 124), 6] <- "Y"
#' RevGNdup[c(111, 115, 128), 7] <- 1
#' RevGNdup[c(113, 117, 120), 7] <- 2
#' RevGNdup[c(116, 119), 7] <- 3
#' RevGNdup[c(123, 125), 7] <- 4
#' RevGNdup[c(126, 127), 7] <- 5
#' subset(RevGNdup, SET_NO==12 & TYPE=="P", select= c(IDKW, DEL, SPLIT))
#' 
#' # Reconstruct ProDup object
#' GNdup2 <- ReconstructProbDup(RevGNdup)
#' lapply(disGNdup, nrow)
#' lapply(GNdup2, nrow)
#' 
#' }
#' @seealso \code{\link[PGRdup]{ProbDup}}, \code{\link[PGRdup]{ReviewProbDup}}
#' @import data.table
#' @importFrom stringi stri_count_fixed
#' @export
ReconstructProbDup <- function(rev) {
  # Check if rev is a data frame
  if (!is.data.frame(rev)) {
    stop("'rev' is not a data frame")
  }
  # Check if core fields are present
  core <- c("SET_NO", "TYPE", "PRIM_ID", "DEL", "SPLIT", "COUNT", "IDKW")
  core2 <- c("K[a]", "K[b]", "K[c]")
  j <- which(core2 %in% colnames(rev))
  core <- union(core, core2[j])
  if (is.element(FALSE, core %in% colnames(rev))) {
    # Check if core fields are present in rev and stop if not
    stop("One or more core fields are missing in 'rev'")
  }
  if (!is.numeric(rev$SET_NO)) {
    stop("'SET_NO' is not of class numeric or integer")
  }
  if (!is.character(rev$TYPE)) {
    stop("'TYPE' is not of class character")
  }
  if (!is.character(rev[, core2[j]])) {
    stop("'K[*]' is not of class character")
  }
  if (!is.character(rev$PRIM_ID)) {
    stop("'PRIM_ID' is not of class character")
  }
  if (!is.character(rev$IDKW)) {
    stop("'IDKW' is not of class character")
  }
  if (!is.numeric(rev$COUNT)) {
    stop("'COUNT' is not of class numeric or integer")
  }
  if (!is.character(rev$DEL)) {
    stop("'DEL' is not of class character")
  }
  if (!is.numeric(rev$SPLIT)) {
    stop("'SPLIT' is not of class numeric or integer")
  }
  # Retrieve method and fields
  method <- regmatches(core2[j], gregexpr("(?<=\\[).*?(?=\\])", core2[j],
                                          perl = TRUE))[[1]]
  fields <- list(k1 = NULL, k2 = NULL)
  fields[[1]] <- colnames(rev)[which(grepl("^K1_", colnames(rev),
                                           perl = TRUE))]
  if (method == "b" | method == "c" ) {
    fields[[2]] <- colnames(rev)[which(grepl("^K2_", colnames(rev),
                                             perl = TRUE))]
  }
  # Remove non-core fields
  #setDT(rev)
  rev <- as.data.table(rev)
  setkey(rev, SET_NO)
  #rev <- data.table(rev, key = "SET_NO")
  rev[ , setdiff(colnames(rev), core) := NULL]
  rev <- rev[ TYPE != ""]
  rev[, SET_NO := as.factor(rev$SET_NO)]
  # Cleanup column DEL
  rev[, DEL := toupper(DEL)]
  rev[, DEL := as.factor(DEL)]
  setattr(rev$DEL,"levels",
          levels(rev$DEL)[levels(rev$DEL) != "Y" & levels(rev$DEL) != "N"] <- "N")
  # Delete sets according to column DEL
  rev <- rev[DEL != "Y"]
  # Cleanup column SPLIT
  rev[, SPLIT := suppressWarnings(as.integer(SPLIT))]
  invisible(rev[is.na(SPLIT), SPLIT := 0])
  rev[, SPLIT := as.factor(SPLIT)]
  # Split sets according to column SPLIT
  rev[, SET_NO := interaction(rev$SET_NO, rev$SPLIT, drop = TRUE)]
  rev[, SET_NO := factor(SET_NO, levels = sort(levels(SET_NO)))]
  setattr(rev$SET_NO,"levels", seq(from = 1, to = length(levels(rev$SET_NO))))
  rev[, SET_NO := as.numeric(SET_NO)]
  setkey(rev, key = "SET_NO")
  # Reconsturct
  rev[, c("DEL", "SPLIT", "COUNT") := NULL]
  rev[, PRIM_ID := paste(get(core2[j]), PRIM_ID, sep = "")]
  rev[, core2[j] := NULL]
  rev <- rev[, list(ID = paste0(sort(unique(PRIM_ID)), collapse = ", "),
                    IDKW = paste0(sort(unique(IDKW)), collapse = ", ")),
             by = c("TYPE", "SET_NO")]
  rev[, COUNT := stri_count_fixed(ID, ",") + 1]
  rev <- rev[COUNT != 1]
  setcolorder(rev, c("SET_NO", "TYPE", "ID", "IDKW", "COUNT"))
  rev <- unique(rev, by = c("TYPE", "ID", "IDKW", "COUNT"))
  out <- list(FuzzyDuplicates = NULL, PhoneticDuplicates = NULL,
              SemanticDuplicates = NULL, DisjointDupicates = NULL)
  attr(out, "method") <- method
  fields[[1]] <-  gsub("^K1_", "", fields[[1]])
  if (method == "b" | method == "c" ) {
    fields[[2]] <- gsub("^K2_", "", fields[[2]])
  }
  attr(out, "fields") <- fields
  rm(fields, method)
  types <- c("F", "P", "S", "D")
  # Reset the SET_NO
  setkey(rev, TYPE, ID, SET_NO)
  rev[, TYPE := as.factor(TYPE)][,SET_NO := as.integer(SET_NO)]
  rev[, SET_NO := 1:.N , by = TYPE]
  rev[, SET_NO := as.numeric(SET_NO)][, TYPE := as.character(TYPE)]
  N <- length(seq_along(out))
  for (i in 1:N) {
    out[[i]] <- setDF(subset(rev, TYPE == types[i]))
    #out[[i]] <- as.data.frame(subset(rev, TYPE == types[i]))
  }
  rev[, TYPE := as.factor(TYPE)]
  if (!"D" %in% levels(rev$TYPE)) {
    out[[4]] <- NULL
  }
  out[which(unlist(lapply(seq_along(out),
                          function(i) dim(out[[i]])[1])) == 0)] <- list(NULL)
  rm(rev)
  class(out) <- "ProbDup"
  return(out)
}
