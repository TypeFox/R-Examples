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



#' Retrieve probable duplicate set information from PGR passport database for 
#' review
#' 
#' \code{ReviewProbDup} retrieves information associated with the probable 
#' duplicate sets from the original PGR passport database(s) from which they 
#' were identified in order to facilitate manual clerical review.
#' 
#' This function helps to retrieve PGR passport information associated with 
#' fuzzy, phonetic or semantic probable duplicate sets in an object of class 
#' \code{ProbDup} from the original databases(s) from which they were 
#' identified. The original information of accessions comprising a set, which 
#' have not been subjected to data standardization can be compared under manual 
#' clerical review for the validation of the set.
#' 
#' By default only the fields(columns) which were used initially for creation of
#' the KWIC indexes using the \code{\link[PGRdup]{KWIC}} function are retrieved.
#' Additional fields(columns) if necessary can be specified using the 
#' \code{extra.db1} and \code{extra.db2} arguments.
#' 
#' The output data frame can be subjected to clerical review either after 
#' exporting into an external spreadsheet using \code{\link[utils]{write.csv}} 
#' function or by using the \code{\link[utils]{edit}} function.
#' 
#' The column \code{DEL} can be used to indicate whether a record has to be 
#' deleted from a set or not. \code{Y} indicates "Yes", and the default \code{N}
#' indicates "No".
#' 
#' The column \code{SPLIT} similarly can be used to indicate whether a record in
#' a set has to be branched into a new set. A set of identical integers in this 
#' column other than the default \code{0} can be used to indicate that they are 
#' to be removed and assembled into a new set.
#' 
#' @param pdup An object of class \code{ProbDup}.
#' @param db1 A data frame of the PGR passport database.
#' @param db2 A data frame of the PGR passport database. Required when 
#'   \code{pdup} was created using more than one KWIC Index.
#' @param extra.db1 A character vector of extra \code{db1} column names to be 
#'   retrieved.
#' @param extra.db2 A character vector of extra \code{db2} column names to be 
#'   retrieved.
#' @param max.count The maximum count of probable duplicate sets whose 
#'   information is to be retrieved.
#' @param insert.blanks logical. If \code{TRUE}, inserts a row of /code{NAs} 
#'   after each set.
#' @return A data frame of the long/narrow form of the probable duplicate sets 
#'   data along with associated fields from the original database(s). The core 
#'   columns in the resulting data frame are as follows: \tabular{ll}{ 
#'   \code{SET_NO} \tab The set number. \cr \code{TYPE} \tab The type of 
#'   probable duplicate set. 'F' for fuzzy, 'P' for phonetic and 'S' for 
#'   semantic matching sets. \cr \code{K[*]} \tab The KWIC index or database of 
#'   origin of the record. The \code{method} is specified within the square 
#'   brackets in the column name.  \cr \code{PRIM_ID} \tab The primary ID of the
#'   accession record from which the set could be identified. \cr \code{IDKW} 
#'   \tab The 'matching' keywords along with the IDs. \cr \code{DEL} \tab Column
#'   to indicate whether record has to be deleted or not. \cr \code{SPLIT} \tab 
#'   Column to indicate whether record has to be branched and assembled into new
#'   set. \cr \code{COUNT} \tab The number of elements in a set. \cr } For the 
#'   retrieved columns(fields) the prefix \code{K*} indicates the KWIC index of 
#'   origin.
#' @note When any primary ID/key records in the fuzzy, phonetic or semantic 
#'   duplicate sets are found to be missing from the original databases 
#'   \code{db1} and \code{db2}, then they are ignored and only the matching 
#'   records are considered for retrieving the information with a warning.
#'   
#'   This may be due to data standardization of the primary ID/key field using 
#'   the function \code{\link[PGRdup]{DataClean}} before creation of the KWIC 
#'   index and subsequent identification of probable duplicate sets. In such a 
#'   case, it is recommended to use an identical data standardization operation 
#'   on the databases \code{db1} and \code{db2} before running this function.
#'   
#'   With \code{R} <= v3.0.2, due to copying of named objects by \code{list()},
#'   \code{Invalid .internal.selfref detected and fixed...} warning can appear,
#'   which may be safely ignored.
#' @seealso \code{\link[PGRdup]{DataClean}}, \code{\link[PGRdup]{KWIC}}, 
#'   \code{\link[PGRdup]{ProbDup}}
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
#' # OR examine and review the duplicate sets after exporting them as a csv file
#' write.csv(file="Duplicate sets for review.csv", x=RevGNdup)
#'                                                                                                                                           
#' }
#' @import data.table
#' @importFrom methods is
#' @export
ReviewProbDup <- function (pdup, db1, db2 = NULL,
                           extra.db1 = NULL, extra.db2 = NULL,
                           max.count = 30, insert.blanks  = TRUE) {
  if (!is(pdup, "ProbDup")) {
    stop('"pdup" is not of class ProbDup')
  }
  method <- attributes(pdup)$method
  fields <- attributes(pdup)$fields
  if (method == "c" | method == "b") {
    if (is.null(db2)) {
      stop(paste("argument 'db2' is missing, with no default.",
                 "\nSecond database is to be specified as method ",
                 method," was used to generate 'pdup'", sep = ""))
    }
  }
  if (is.element(FALSE, fields[[1]] %in% colnames(db1))) {
    # Check if fields are present in db1 and stop if not
    stop("One or more kwic1 fields are missing in 'db1'")
  }
  if (!is.null(extra.db1) && !is.vector(extra.db1, mode = "character")) {
    stop('"extra.db1" is not a character vector')
  }
  if (is.element(FALSE, extra.db1 %in% colnames(db1))) {
    # Check if extra fields are present in db1 and stop if not
    warning(paste("One or more extra kwic1 fields are missing in 'db1'.",
                  "Only the default kwic1 fields will be retrieved"))
  }
  fields[[1]] <- union(fields[[1]], extra.db1)
  #setDT(db1)
  db1 <- as.data.table(db1)
  if (!identical(setdiff(colnames(db1), fields[[1]]), character(0))){
    db1[, setdiff(colnames(db1), fields[[1]]) := NULL]
  }
  db1[, K1_PRIM_ID := get(fields[[1]][1])]
  setcolorder(db1, neworder = union("K1_PRIM_ID",
                                    setdiff(colnames(db1), "K1_PRIM_ID")))
  db1[, K1_PRIM_ID := as.character(K1_PRIM_ID)]
  setkey(db1, "K1_PRIM_ID")
  if (method == "c" | method == "b") {
    if (is.element(FALSE, fields[[2]] %in% colnames(db2)) == TRUE) {
      # Check if fields are present in db2 and stop if not
      stop("One or more kwic2 fields are missing in 'db2'")
    }
    if (!is.null(extra.db2) && !is.vector(extra.db2, mode = "character")) {
      stop("'extra.db2' is not a character vector")
    }
    if (is.element(FALSE, extra.db2 %in% colnames(db2)) == TRUE) {
      # Check if extra fields are present in db2 and stop if not
      warning(paste("One or more extra kwic2 fields are missing in 'db2'.",
                    "Only the default kwic2 fields will be retrieved"))
    }
    fields[[2]] <- union(fields[[2]], extra.db2)
    #setDT(db2)
    db2 <- as.data.table(db2)
    if (!identical(setdiff(colnames(db2), fields[[1]]), character(0))) {
      db2[, setdiff(colnames(db2), fields[[2]]) := NULL]
    }
    db2[, K2_PRIM_ID := get(fields[[2]][1])]
    setcolorder(db2, neworder = union("K2_PRIM_ID",
                                      setdiff(colnames(db2), "K2_PRIM_ID")))
    db2[, K2_PRIM_ID := as.character(K2_PRIM_ID)]
    setkey(db2, "K2_PRIM_ID")
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
                                                         PRIM_ID = gsub(":.*", "", V1))]
      # Create K1_PRIM_ID column for merging with db1
      pdup[[i]][, K1_PRIM_ID := gsub("\\[K2\\].*", "", PRIM_ID, perl = TRUE)]
      pdup[[i]][, K1_PRIM_ID := gsub("\\[K1\\]", "", K1_PRIM_ID, perl = TRUE)]
      pdup[[i]][, K1_PRIM_ID := as.character(K1_PRIM_ID)]
      setkey(pdup[[i]],K1_PRIM_ID)
      # Check for IDs in pdup missing in db1
      if (length(setdiff(pdup[[i]]$K1_PRIM_ID[pdup[[i]]$K1_PRIM_ID != ""],
                         db1$K1_PRIM_ID)) != 0) {
        warning(paste("Encountered primary ID records in",
                      types2[i], "probable duplicate sets missing from 'db1'.",
                      "\nOnly matching records are merged"))
      }
      # Merge with db1
      pdup[[i]] <- merge(pdup[[i]],db1, all.x = TRUE)
      pdup[[i]][, K1_PRIM_ID := NULL]
      setkey(pdup[[i]],"SET_NO")
      # Set prefix to merged fields
      setnames(pdup[[i]], old = fields[[1]],
               new =  paste("K1", fields[[1]], sep = "_"))
      if (method == "c" | method == "b") {
        # Create K2_PRIM_ID column for merging with db2
        pdup[[i]][, K2_PRIM_ID := gsub("\\[K1\\].*", "", PRIM_ID, perl = TRUE)]
        pdup[[i]][, K2_PRIM_ID := gsub("\\[K2\\]", "", K2_PRIM_ID, perl = TRUE)]
        pdup[[i]][, K2_PRIM_ID := as.character(K2_PRIM_ID)]
        setkey(pdup[[i]],K2_PRIM_ID)
        # Check for IDs in pdup missing in db1
        if (length(setdiff(pdup[[i]]$K2_PRIM_ID[pdup[[i]]$K2_PRIM_ID != ""],
                           db2$K2_PRIM_ID)) != 0) {
          warning(paste("Encountered primary ID records in",
                        types2[i],
                        "probable duplicate sets missing from 'db2'.",
                        "\nOnly matching records are merged"))
        }
        # Merge with db1
        pdup[[i]] <- merge(pdup[[i]],db2, all.x = TRUE)
        pdup[[i]][, K2_PRIM_ID := NULL]
        setkey(pdup[[i]],"SET_NO")
        setnames(pdup[[i]], old = fields[[2]],
                 new =  paste("K2", fields[[2]], sep = "_"))
      }
    }
  }
  # rbind pdup list
  pdup <- rbindlist(pdup)
  # Split K* from PRIM_ID column
  pdup[, PRIM_ID := gsub("\\]", "\\]_", PRIM_ID, perl = TRUE)]
  pdup[, K := gsub("_.*", "", PRIM_ID, perl = TRUE)]
  pdup[, PRIM_ID := gsub("\\[K1\\]_", "", PRIM_ID, perl = TRUE)]
  pdup[, PRIM_ID := gsub("\\[K2\\]_", "", PRIM_ID, perl = TRUE)]
  # Add review columns
  pdup[,DEL := "N" ]
  pdup[,SPLIT := 0 ]
  # Reset column order
  nameslist <- union(c("SET_NO", "TYPE", "K", "PRIM_ID",
                       "IDKW", "DEL", "SPLIT", "COUNT"),
                    colnames(pdup))
  setcolorder(x = pdup, neworder = nameslist)
  # Add prefix to additional fields
  if (!is.null(extra.db1) && !is.element(FALSE,
                                        paste("K1", extra.db1, sep = "_") %in% colnames(pdup))) {
    setnames(pdup, old = paste("K1", extra.db1, sep = "_"),
             new = paste("K1X", extra.db1, sep = "_"))
  }
  if (method == "c" | method == "b") {
    if (!is.null(extra.db2) && !is.element(FALSE,
                                          paste("K2", extra.db2,
                                                sep = "_") %in% colnames(pdup))) {
      setnames(pdup, old = paste("K2", extra.db2, sep = "_"),
               new = paste("K2X", extra.db2, sep = "_"))
    }
  }
  # Insert blanks
  setkey(pdup, SET_NO)
  setkey(pdup, NULL)
  if (insert.blanks  == TRUE) {
    pdup[, TEMP := as.factor(SET_NO)]
    pdup[, TEMP := interaction(pdup$TEMP, as.factor(pdup$TYPE), drop = TRUE)]
    setattr(pdup$TEMP,"levels", seq(from = 1, to = length(levels(pdup$TEMP))))
    pdup[, TEMP := as.numeric(TEMP)]
    pdup <- setDT(pdup)[pdup[, c(.I, NA), TEMP]$V1][!.N]
    pdup[, TEMP := NULL]
    pdup[is.na(TYPE), TYPE := ""]
  }
  setnames(pdup, old = "K", new = paste("K[", method, "]", sep = ""))
  setDF(pdup)
  return(pdup)
}
