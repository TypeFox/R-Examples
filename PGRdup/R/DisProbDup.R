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



#' Get disjoint probable duplicate sets
#' 
#' \code{DisProbDup} finds and joins intersecting sets in an object of class 
#' \code{ProbDup} to get disjoint probable duplicate sets.
#' 
#' This function considers the accession primary keys/IDs for finding 
#' intersecting sets and subsequently joins them to retrieve disjoint sets. 
#' These operations are implemented utilizing the \code{\link[igraph]{igraph}} 
#' package functions.
#' 
#' Disjoint sets are retrieved either individually for each type of probable 
#' duplicate sets or considering all type of sets simultaneously. In case of the
#' latter, the disjoint of all the type of sets alone are returned in the output
#' as an additional data frame \code{DisjointDuplicates} in an object of class 
#' \code{ProbDup}
#' 
#' @param pdup An object of class \code{ProbDup}.
#' @param combine A character vector indicating the type of sets to be 
#'   considered together for retrieving disjoint sets. If \code{NULL}, then 
#'   disjoint sets within each type are retrieved (see \strong{Details}).
#' @return Returns an object of class \code{ProbDup} with either the disjoint 
#'   sets within each type - \code{FuzzyDuplicates}, \code{PhoneticDuplicates} 
#'   and \code{SemanticDuplicates} when \code{combine = NULL} or the combined 
#'   disjoint duplicate sets as an additional element \code{DisjointDupicates}
#'   according to the choice specified in the argument \code{combine}.
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
#' lapply(GNdup, dim)
#' 
#' # Get disjoint probable duplicate sets of each kind
#' disGNdup1 <- DisProbDup(GNdup, combine = NULL)
#' lapply(disGNdup1, nrow)
#' 
#' # Get disjoint probable duplicate sets combining all the kinds of sets
#' disGNdup2 <- DisProbDup(GNdup, combine = c("F", "P", "S"))
#' lapply(disGNdup2, nrow)
#'                   
#' }                 
#' @seealso \code{\link[PGRdup]{ProbDup}}
#' @import igraph
#' @import data.table
#' @importFrom methods is
#' @importFrom stats embed
#' @importFrom utils stack
#' @importFrom stringi stri_count_fixed
#' @export
DisProbDup <- function(pdup, combine = c("F", "P", "S")) {
  # Check arguments
  if (!is(pdup, "ProbDup")) {
    stop('"pdup" is not of class ProbDup')
  }
  out <- list(FuzzyDuplicates = NULL, PhoneticDuplicates = NULL,
              SemanticDuplicates = NULL, DisjointDupicates = NULL)
  attr(out, "method") <- attr(pdup, "method")
  attr(out, "fields") <- attr(pdup, "fields")
  types <- c("F", "P", "S", "D")
  # Check if Disjoint duplicates is in list and convert combine and pdup
  if (!is.null(pdup[["DisjointDupicates"]])) {
    pdup[1:3] <- NULL
    combine <- NULL
    warning(paste("Only Disjoint duplicate sets encountered in 'pdup'",
                  "\nFurther disjoint sets, if any are returned", sep = ""))
  }
  cond <- length(!unlist(lapply(pdup[1:3],
                                is.null))[!unlist(lapply(pdup[1:3],
                                                         is.null))]) == 1
  if (is.null(pdup[["DisjointDupicates"]]) & cond) {
    combine <- NULL
    warning(paste("Only one kind of probable duplicate set encountered in 'pdup'",
                  "\nDisjoint sets within the same are returned", sep = ""))
  }
  if (!is.null(combine)) {
    combine <- match.arg(combine, c("F", "P", "S"), several.ok = TRUE)
    if (length(combine == 1) & is.element(TRUE,
                                          lapply(pdup[which(types %in% combine)],
                                                 is.null))) {
      p <- unlist(lapply(pdup[which(types %in% combine)], is.null))
      stop(paste("The following set specified in argument 'combine' is missing from 'pdup'",
                 "\n", paste("#",names(p[p]), collapse = "\n")))
    }
    # Check and report any missing sets specified in combine
    if (is.element(TRUE, lapply(pdup[which(types %in% combine)], is.null))) {
      p <- unlist(lapply(pdup[which(types %in% combine)], is.null))
      warning(paste("The following set specified in argument 'combine' is missing from 'pdup'",
                    "\n", paste("#",names(p[p]), collapse = "\n"),
                    "\nJoint disjoint of only the remaining sets are returned",
                    sep = ""))
    }
    # Rbind all in combine into pdup[4] and assign NULL to rest
    pdup[which(!types %in% combine)] <- list(NULL)
    pdup[[4]] <- rbindlist(pdup[1:3])
    pdup[which(types %in% combine)] <- list(NULL)
    pdup[[4]][, TYPE := "D"]
  }
  N <- length(seq_along(pdup))
  for (i in 1:N) {
    if (!is.null(pdup[[i]])) {
      # Get and cast the disjoint sets by PRIM_ID
      idlist <- strsplit(pdup[[i]]$ID, ", ")
      idcomb <- do.call("rbind",lapply(idlist, embed, 2))
      gg <- graph.edgelist(idcomb, directed = FALSE)
      disidlist <- split(V(gg)$name, clusters(gg)$membership)
      disidlist <- stack(lapply(disidlist, sort))
      #setDT(disidlist)
      disidlist <- as.data.table(disidlist)
      setkey(disidlist, ind)
      #disidlist <- data.table(stack(disidlist), key = "ind")
      disidlist[, ind := as.integer(ind)]
      setnames(disidlist, old = c("values", "ind"), new = c("ID", "SET_NO"))
      setkey(disidlist, ID)
      rm(gg, idlist, idcomb)
      #setDT(pdup[[i]])
      pdup[[i]] <- as.data.table(pdup[[i]])
      #pdup[[i]] <- data.table(pdup[[i]])
      pdup[[i]] <- pdup[[i]][, .(unlist(strsplit(IDKW, ", ", TRUE))), by = TYPE][,
                               .(IDKW = toString(sort(unique(unlist(strsplit(V1,", ")))))), .(TYPE, ID = gsub(":.*", "", V1))]
      pdup[[i]] <- unique(pdup[[i]])
      setkey(pdup[[i]], ID)
      disidlist <- merge(disidlist, pdup[[i]])
      disidlist <- disidlist[, list(ID = paste0(sort(unique(ID)),
                                                collapse = ", "),
                                    IDKW = paste0(sort(unique(IDKW)),
                                                  collapse = ", ")),
                             by = c("TYPE", "SET_NO")]
      disidlist[, COUNT := stri_count_fixed(ID, ",") + 1]
      setkey(disidlist, ID)
      disidlist[, SET_NO := seq(1, nrow(disidlist))]
      disidlist[, SET_NO := as.numeric(SET_NO)]
      setkey(disidlist, SET_NO)
      setcolorder(disidlist, c("SET_NO", "TYPE", "ID", "IDKW", "COUNT"))
      out[[i]] <- setDF(disidlist)
      #out[[i]] <- as.data.frame(disidlist)
    }
  }
  if (is.null(combine)) {
    out[[4]] <- NULL
  }
  class(out) <- "ProbDup"
  return(out)
}
