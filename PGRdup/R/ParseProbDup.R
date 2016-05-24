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



#' Parse an object of class \code{ProbDup} to a data frame.
#' 
#' \code{ParseProbDup} converts an object of class \code{ProbDup} to a data frame for export.
#' 
#' @param pdup An object of class \code{ProbDup}.
#' @param max.count The maximum count of probable duplicate sets which are to be parsed to a data frame.
#' @param insert.blanks logical. If \code{TRUE}, inserts a row of \code{NAs} 
#'   after each set.
#' @return A data frame of the long/narrow form of the probable duplicate sets 
#'   data with the following core columns: \tabular{ll}{ 
#'   \code{SET_NO} \tab The set number. \cr \code{TYPE} \tab The type of 
#'   probable duplicate set. 'F' for fuzzy, 'P' for phonetic and 'S' for 
#'   semantic matching sets. \cr \code{K} \tab The KWIC index or database of 
#'   origin of the record. The \code{method} is specified within the square 
#'   brackets in the column name.  \cr \code{PRIM_ID} \tab The primary ID of the
#'   accession record from which the set could be identified. \cr \code{IDKW} 
#'   \tab The 'matching' keywords along with the IDs. \cr \code{COUNT} \tab The number of elements in a set. \cr } For the 
#'   retrieved columns(fields) the prefix \code{K*} indicates the KWIC index of 
#'   origin.
#' @examples
#' \dontrun{
#' 
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
#' # Convert to data frame of sets               
#' GNdupParsed <- ParseProbDup(GNdup)
#' 
#' }
#' @seealso 
#'   \code{\link[PGRdup]{ProbDup}},
#' @import data.table
#' @importFrom methods is
#' @export
ParseProbDup <- function(pdup, max.count = 30,
                          insert.blanks  = TRUE) {
  if (!is(pdup, "ProbDup")) {
    stop('"pdup" is not of class ProbDup')
  }
  N <- length(seq_along(pdup))
  types <- c("F", "P", "S", "D")
  types2 <- c("Fuzzy", "Phonetic", "Semantic", "Disjoint")
  for (i in 1:N) {
    if (!is.null(pdup[[i]])) {
      #setDT(pdup[[i]])
      pdup[[i]] <- as.data.table(pdup[[i]])
      pdup[[i]] <- subset(pdup[[i]], COUNT <= max.count)
      pdup[[i]] <- unique(pdup[[i]])
      # Reset SET_NO to take into account deleted sets with coutn > max.count
      pdup[[i]][, Seq := 1:.N]
      # Cast ID and IDKW by SET_NO
      pdup[[i]] <- pdup[[i]][, .(unlist(strsplit(IDKW, ", ", TRUE))),
                             by = list(SET_NO, TYPE, COUNT)][,
                               .(IDKW = toString(V1)), .(SET_NO, TYPE, COUNT,
                                                         PRIM_ID = gsub(":.*",
                                                                        "",
                                                                        V1))]
    }
  }
  # rbind pdup list
  pdup <- rbindlist(pdup)
  # Split K* from PRIM_ID column
  pdup[,PRIM_ID := gsub("\\]", "\\]_", PRIM_ID, perl = TRUE)]
  pdup[,K := gsub("_.*", "", PRIM_ID, perl = TRUE)]
  pdup[,PRIM_ID := gsub("\\[K1\\]_", "", PRIM_ID, perl = TRUE)]
  pdup[,PRIM_ID := gsub("\\[K2\\]_", "", PRIM_ID, perl = TRUE)]
  # Reset column order
  nameslist <- c("SET_NO", "TYPE", "K", "PRIM_ID", "IDKW", "COUNT")
  setcolorder(x = pdup, neworder = nameslist)
  # Insert blanks
  setkey(pdup,NULL)
  if (insert.blanks  == TRUE) {
    pdup[,TEMP := as.factor(SET_NO)]
    pdup[,TEMP := interaction(pdup$TEMP, as.factor(pdup$TYPE), drop = TRUE)]
    setattr(pdup$TEMP,"levels", seq(from = 1, to = length(levels(pdup$TEMP))))
    pdup[,TEMP := as.numeric(TEMP)]
    pdup <- setDT(pdup)[pdup[, c(.I, NA), TEMP]$V1][!.N]
    pdup[, TEMP := NULL]
    pdup[is.na(TYPE), TYPE := ""]
  }
  setDF(pdup)
  return(pdup)
}
