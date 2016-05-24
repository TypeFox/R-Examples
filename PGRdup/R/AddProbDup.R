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



#' Add probable duplicate sets fields to the PGR passport database
#' 
#' \code{AddProbDup} adds the fuzzy, phonetic and semantic probable duplicates 
#' sets data fields from an object of class \code{ProbDup} to the original PGR 
#' passport database.
#' 
#' This function helps to add information associated with identified fuzzy, 
#' phonetic and semantic probable duplicate sets using the 
#' \code{\link[PGRdup]{ProbDup}} function to the original PGR passport database.
#' Associated data fields such as \code{SET_NO}, \code{ID} and \code{IDKW} are 
#' added based on the \code{PRIM_ID} field(column).
#' 
#' In case more than one KWIC index was used to generate the object of class 
#' \code{ProbDup}, the argument \code{addto} can be used to specify to which 
#' database the data fields are to be added. The default \code{"I"} indicates 
#' the database from which the first KWIC index was created and \code{"II"} 
#' indicates the database from which the second index was created.
#' 
#' @param pdup An object of class \code{ProbDup}.
#' @param db A data frame of the PGR passport database.
#' @param addto Either \code{"I"} or \code{"II"} indicating the database to 
#'   which the data.fields are to be added (see \strong{Details}).
#' @param max.count The maximum count of probable duplicate sets whose 
#'   information is to be retrieved.
#' @return A data frame of the PGR passport database with the probable duplicate
#'   sets fields added.
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
#' # Add the duplicates sets to the original database                 
#' GNwithdup <-  AddProbDup(pdup = GNdup, db = GN1000, addto = "I")                  
#' 
#' }
#' @seealso \code{\link[PGRdup]{DataClean}}, \code{\link[PGRdup]{KWIC}},
#'   \code{\link[PGRdup]{ProbDup}}
#' @note When any primary ID/key records in the fuzzy, phonetic or semantic 
#'   duplicate sets are found to be missing from the original database 
#'   \code{db}, then they are ignored and only the matching records are 
#'   considered for adding the information with a warning.
#'   
#'   This may be due to data standardization of the primary ID/key field using 
#'   the function \code{\link[PGRdup]{DataClean}} before creation of the KWIC 
#'   index and subsequent identification of probable duplicate sets. In such a 
#'   case, it is recommended to use an identical data standardization operation 
#'   on the database \code{db} before running this function.
#' @import data.table
#' @importFrom methods is
#' @export
AddProbDup <- function(pdup, db, addto = c("I", "II"), max.count = 30) {
  if (!is(pdup, "ProbDup")) {
    stop('"pdup" is not of class ProbDup')
  }
  addto <- match.arg(addto, several.ok = FALSE)
  method <- attributes(pdup)$method
  if (addto == "I") {
    p <- "[K1]"
    fields <- attributes(pdup)$fields[[1]]
  }
  if (addto == "II") {
    if (is.null(attributes(pdup)$fields[[2]])) {
      stop("argument 'addto = 'II', but pdup was not generated using a second KWIC index" )
    }
    p <- "[K2]"
    fields <- attributes(pdup)$fields[[2]]
  }
  if (is.data.frame(db) == FALSE) {
    # Check if db is a data frame and stop if not
    stop('"db" is not a data frame')
  }
  if (!fields[1] %in% colnames(db)) {
    stop(paste('Column(field) matching the primary ID/key "',
               fields[1],'" is not found in "db"', sep = ""))
  }
  if (is.element("FALSE", fields %in% colnames(db)) == TRUE) {
    # Check if fields are present in x and stop if not
    warning('One or more fields used for generating "pdup" are missing in "db"')
  }
  #setDT(db)
  db <- as.data.table(db)
  setkeyv(db, fields[[1]])
  # Retrieve the sets with "[K*]" according to argument addto
  pdup <- lapply(pdup, function(x) x[grepl(p , x$ID, fixed = TRUE),])
  #invisible(lapply(seq_along(pdup), function(i) pdup[[i]]$PRIM_ID <<- gsub(p, "", pdup[[i]]$PRIM_ID, fixed = TRUE)))
  types <- c("F", "P", "S", "D")
  types2 <- c("Fuzzy", "Phonetic", "Semantic", "Disjoint")
  N <- length(pdup)
  for (i in 1:N) {
    if (!is.null(pdup[[i]])) {
      #setDT(pdup[[i]])
      pdup[[i]] <- as.data.table(pdup[[i]])
      pdup[[i]] <- subset(pdup[[i]], COUNT <= max.count)
      pdup[[i]][, TYPE := NULL]
      pdup[[i]] <- pdup[[i]][, list(PRIM_ID = unlist(strsplit(ID , ", " ))) ,
                             by = list(SET_NO, IDKW, ID)]
      if (method == "b" | method == "c") {
        if (addto == "I") {
          pdup[[i]] <- pdup[[i]][grepl("\\[K1\\].*", PRIM_ID, perl = TRUE)]
        }
        if (addto == "II") {
          pdup[[i]] <- pdup[[i]][grepl("\\[K2\\].*", PRIM_ID, perl = TRUE)]
        }
      }
      pdup[[i]][, PRIM_ID := gsub(p, "", pdup[[i]]$PRIM_ID, fixed = TRUE)]
      setkey(pdup[[i]], PRIM_ID)
      # Aggregate according to ID
      pdup[[i]] <- pdup[[i]][, list(SET_NO = paste0(sort(unique(SET_NO)),
                                                    collapse = ", "),
                                    IDKW = paste0(sort(unique(unlist(strsplit(IDKW , ", ")))),
                                                  collapse = ", "),
                                    ID = paste0(sort(unique(unlist(strsplit(ID, ", ")))),
                                                collapse = ", ")), by = PRIM_ID]
      setkey(pdup[[i]], PRIM_ID)
      setcolorder(pdup[[i]], c("PRIM_ID", "SET_NO", "ID", "IDKW"))
      setnames(pdup[[i]], old = colnames(pdup[[i]]),
               new = paste(types[i], colnames(pdup[[i]]), sep = "_"))
      setnames(pdup[[i]], old = paste(types[i], "PRIM_ID", sep = "_"),
               new = fields[1])
      if (length(setdiff(pdup[[i]][,get(fields[1])],
                         db[,get(fields[1])])) != 0) {
        warning(paste("Encountered primary ID records in",
                      types2[i], "probable duplicate sets missing from 'db'.",
                      "\nOnly matching records are merged"))
      }
      db <- merge(db, pdup[[i]], all.x = TRUE)
    }
  }
  rm(pdup)
  db[is.na(db)]   <- ""
  setDF(db)
  return(db)
}
