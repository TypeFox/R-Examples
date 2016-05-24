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



#' Identify Probable Duplicates of Accessions
#' 
#' \code{ProbDup} identifies probable duplicates of germplasm accessions in KWIC
#' indexes created from PGR passport databases using fuzzy, phonetic and 
#' semantic matching strategies.
#' 
#' This function performs fuzzy, phonetic and semantic matching of keywords in 
#' KWIC indexes of PGR passport databases (created using 
#' \code{\link[PGRdup]{KWIC}} function) to identify probable duplicates of 
#' germplasm accessions. The function can execute matching according to either 
#' of the following three methods as specified by the \code{method} argument.
#' 
#' \describe{ \item{Method \code{a}:}{Perform string matching of keywords in a 
#' single KWIC index to identify probable duplicates of accessions in a single 
#' PGR passport database.} \item{Method \code{b}:}{Perform string matching of 
#' keywords in the first KWIC index (query) with that of the keywords in the 
#' second index (source) to identify probable duplicates of accessions of the 
#' first PGR passport database among the accessions in the second database.} 
#' \item{Method \code{c}:}{Perform string matching of keywords in two different 
#' KWIC indexes jointly to identify probable duplicates of accessions from among
#' two PGR passport databases.}}
#' 
#' \strong{Fuzzy matching} or approximate string matching of keywords is carried
#' out by computing the generalized levenshtein (edit) distance between them. 
#' This distance measure  counts the number of deletions, insertions and 
#' substitutions necessary to turn one string to the another. A distance of up 
#' to \code{max.dist} are considered for a match.
#' 
#' Exact matching will be enforced when the argument \code{force.exact} is 
#' \code{TRUE}. It can be used to avoid fuzzy matching when the number of 
#' alphabet characters in keywords is lesser than a critical value 
#' (\code{max.alpha}). Similarly, the value of \code{max.digit} can also be set 
#' according to the requirements. The default value of \code{Inf} avoids fuzzy 
#' matching and enforces exact matching for all keywords having any numerical 
#' characters. If \code{max.digit} and \code{max.alpha} are both set to 
#' \code{Inf}, exact matching will be enforced for all the keywords.
#' 
#' When exact matching is enforced, for keywords having both alphabet and 
#' numeric characters and with the number of alphabet characters greater than 
#' \code{max.digit}, matching will be carried out separately for alphabet and 
#' numeric characters present.
#' 
#' \strong{Phonetic matching} of keywords is carried out using the Double 
#' Metaphone phonetic algorithm (\code{\link[PGRdup]{DoubleMetaphone}}) to 
#' identify keywords that have the similar pronunciation. Either the 
#' \code{primary} or \code{alternate} encodings can be used by specifying the 
#' \code{encoding} argument. The argument \code{phon.min.alpha} sets the limits 
#' for the number of alphabet characters to be present in a string for executing
#' phonetic matching. Similarly \code{min.enc} sets the limits for the number of
#' characters to be present in the encoding of a keyword for phonetic matching.
#' 
#' \strong{Semantic matching} matches keywords based on a list of accession name
#' synonyms supplied as list with character vectors of synonym sets (synsets) to
#' the \code{syn} argument. Synonyms in this context refers to interchangeable 
#' identifiers or names by which an accession is recognized. Multiple keywords
#' specified as members of the same synset in \code{syn} are merged together.
#' To facilitate accurate identification of synonyms from the KWIC index,
#' identical data standardization operations using the 
#' \code{\link[PGRdup]{MergeKW}} and \code{\link[PGRdup]{DataClean}} functions 
#' for both the original database fields and the synset list are recommended.
#' 
#' The probable duplicate sets identified initially here may be intersecting 
#' with other sets. To get the disjoint sets after the union of all the 
#' intersecting sets use the \code{\link[PGRdup]{DisProbDup}} function.
#' 
#' The function \code{\link[PGRdup]{AddProbDup}} can be used to add the 
#' information associated with the identified sets in an object of class 
#' \code{ProbDup} as fields(columns) to the original PGR passport database.
#' 
#' All of the string matching operations here are executed through the 
#' \code{\link[stringdist]{stringdist-package}} functions.
#' 
#' @param kwic1 An object of class \code{KWIC}.
#' @param kwic2 An object of class \code{KWIC}. Required for \code{method} 
#'   \code{"b"} and \code{"c"} only (see \strong{Details}).
#' @param method The method to be followed for identification of probable 
#'   duplicates. Either \code{"a"}, \code{"b"} or \code{"c"}. (see 
#'   \strong{Details}).
#' @param chunksize A value indicating the size of KWIC index keyword block to 
#'   be used for searching for matches at a time in case of large number of 
#'   keywords(see \strong{Note}).
#' @param useBytes logical. If \code{TRUE}, performs byte-wise comparison 
#'   instead of character-wise comparison (see \strong{Note}).
#' @param excep A vector of the keywords in KWIC not to be used for probable 
#'   duplicate search (see \strong{Details}).
#' @param fuzzy logical. If \code{TRUE} identifies probable duplicates based on 
#'   fuzzy matching.
#' @param max.dist The maximum levenshtein distance between keyword strings 
#'   allowed for a match. Default is 3 (see \strong{Details}).
#' @param force.exact logical. If \code{TRUE}, enforces exact matching instead 
#'   of fuzzy matching for keyword strings which match the criteria specified in
#'   arguments \code{max.alpha} and \code{max.digit} (see \strong{Details}).
#' @param max.alpha Maximum number of alphabet characters present in a keyword 
#'   string up to which exact matching is enforced rather than fuzzy matching. 
#'   Default is 4 (see \strong{Details}).
#' @param max.digit Maximum number of numeric characters present in a keyword 
#'   string up to which exact matching is enforced rather than fuzzy matching. 
#'   Default is Inf (see \strong{Details}).
#' @param phonetic logical. If \code{TRUE} identifies probable duplicates based 
#'   on phonetic matching.
#' @param encoding Double metaphone encoding for phonetic matching. The default 
#'   is \code{"primary"} (see \strong{Details}).
#' @param phon.min.alpha Minimum number of alphabet characters to be present in 
#'   a keyword string for phonetic matching (see \strong{Details}).
#' @param min.enc Minimum number of characters to be be present in double 
#'   metaphone encoding of a keyword string for phonetic matching (see 
#'   \strong{Details}).
#' @param semantic logical. If \code{TRUE} identifies probable duplicates based 
#'   on semantic matching.
#' @param syn A list with character vectors of synsets (see \strong{Details}).
#' @note As the number of keywords in the KWIC indexes increases, the memory 
#'   consumption by the function also increases. For string matching, this 
#'   function relies upon creation of a \eqn{n}*\eqn{m} matrix of all possible 
#'   keyword pairs for comparison, where \eqn{n} and \eqn{m} are the number of 
#'   keywords in the query and source indexes respectively. This can lead to 
#'   \code{cannot allocate vector of size} errors in case very large KWIC 
#'   indexes where the comparison matrix is too large to reside in memory. In 
#'   such a case, try to adjust the \code{chunksize} argument to get the 
#'   appropriate size of the KWIC index keyword block to be used for searching 
#'   for matches at a time. However a smaller chunksize may lead to longer 
#'   computation time due to the memory-time trade-off.
#'   
#'   The progress of matching is displayed in the console as number of blocks 
#'   completed out of total (e.g. 6 / 30), the percentage of achievement (e.g. 
#'   30\%) and a text-based progress bar.
#'   
#'   In case of multi-byte characters in keywords, the matching speed is further
#'   dependent upon the \code{useBytes} argument as described in 
#'   \strong{Encoding issues} for the \code{\link[stringdist]{stringdist}} 
#'   function, which is made use of here for string matching.
#' @references van der Loo, M. P. J. (2014). The \code{stringdist} Package for 
#'   Approximate String Matching. R Journal, 6(1), 111-122.
#' @examples
#' \dontrun{
#' 
#' # Method "a"
#' #===========
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
#' GNdup
#' 
#' # Method "b and c"
#' #=================
#' 
#' # Load PGR passport databases
#' GN1 <- GN1000[!grepl("^ICG", GN1000$DonorID), ]
#' GN1$DonorID <- NULL
#' GN2 <- GN1000[grepl("^ICG", GN1000$DonorID), ]
#' GN2 <- GN2[!grepl("S", GN2$DonorID), ]
#' GN2$NationalID <- NULL
#' 
#' # Specify as a vector the database fields to be used
#' GN1fields <- c("NationalID", "CollNo", "OtherID1", "OtherID2")
#' GN2fields <- c("DonorID", "CollNo", "OtherID1", "OtherID2")
#' 
#' # Clean the data
#' GN1[GN1fields] <- lapply(GN1[GN1fields], function(x) DataClean(x))
#' GN2[GN2fields] <- lapply(GN2[GN2fields], function(x) DataClean(x))
#' y1 <- list(c("Gujarat", "Dwarf"), c("Castle", "Cary"), c("Small", "Japan"),
#' c("Big", "Japan"), c("Mani", "Blanco"), c("Uganda", "Erect"),
#' c("Mota", "Company"))
#' y2 <- c("Dark", "Light", "Small", "Improved", "Punjab", "SAM")
#' y3 <- c("Local", "Bold", "Cary", "Mutant", "Runner", "Giant", "No.",
#'         "Bunch", "Peanut")
#' GN1[GN1fields] <- lapply(GN1[GN1fields], function(x) MergeKW(x, y1, delim = c("space", "dash")))
#' GN1[GN1fields] <- lapply(GN1[GN1fields], function(x) MergePrefix(x, y2, delim = c("space", "dash")))
#' GN1[GN1fields] <- lapply(GN1[GN1fields], function(x) MergeSuffix(x, y3, delim = c("space", "dash")))
#' GN2[GN2fields] <- lapply(GN2[GN2fields], function(x) MergeKW(x, y1, delim = c("space", "dash")))
#' GN2[GN2fields] <- lapply(GN2[GN2fields], function(x) MergePrefix(x, y2, delim = c("space", "dash")))
#' GN2[GN2fields] <- lapply(GN2[GN2fields], function(x) MergeSuffix(x, y3, delim = c("space", "dash")))
#' 
#' # Remove duplicated DonorID records in GN2
#' GN2 <- GN2[!duplicated(GN2$DonorID), ]
#' 
#' # Generate KWIC index
#' GN1KWIC <- KWIC(GN1, GN1fields)
#' GN2KWIC <- KWIC(GN2, GN2fields)
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
#' GNdupb <- ProbDup(kwic1 = GN1KWIC, kwic2 = GN2KWIC, method = "b",
#'                   excep = exep, fuzzy = TRUE, phonetic = TRUE,
#'                   encoding = "primary", semantic = TRUE, syn = syn)
#' GNdupb
#'                   
#' GNdupc <- ProbDup(kwic1 = GN1KWIC, kwic2 = GN2KWIC, method = "c",
#'                   excep = exep, fuzzy = TRUE, phonetic = TRUE,
#'                   encoding = "primary", semantic = TRUE, syn = syn)
#' GNdupc
#' 
#' }
#' @seealso \code{\link[PGRdup]{KWIC}}, \code{\link[PGRdup]{DoubleMetaphone}} 
#'   \code{\link[stringdist]{stringdistmatrix}}, \code{\link[utils]{adist}}, 
#'   \code{\link[PGRdup]{print.ProbDup}}
#' @return A list of class \code{ProbDup} containing the following data frames 
#'   of probable duplicate sets identified along with the corresponding keywords
#'   and set counts: \enumerate{ \item \code{FuzzyDuplicates} \item 
#'   \code{PhoneticDuplicates} \item \code{SemanticDuplicates} } Each data frame
#'   has the following columns: \tabular{ll}{ \code{SET_NO} \tab The set number.
#'   \cr \code{TYPE} \tab The type of probable duplicate set. 'F' for fuzzy, 'P'
#'   for phonetic and 'S' for semantic matching sets. \cr \code{ID} \tab The 
#'   primary IDs of records of accessions comprising a set. \cr \code{ID:KW} 
#'   \tab The 'matching' keywords along with the IDs. \cr \code{COUNT} \tab The 
#'   number of elements in a set. \cr }
#'   
#'   The prefix \code{[K*]} indicates the KWIC index of origin of the KEYWORD or
#'   PRIM_ID.
#' @import igraph
#' @import stringdist
#' @import data.table
#' @importFrom stringi stri_count_fixed
#' @importFrom stringi stri_count_regex
#' @importFrom methods is
#' @importFrom utils stack
#' @importFrom utils capture.output
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom stats embed
#' @export ProbDup
#' @export print.ProbDup
ProbDup <- function (kwic1, kwic2 = NULL, method = c("a", "b", "c"),
                     excep = NULL, chunksize = 1000, useBytes = TRUE,
                     fuzzy = TRUE, max.dist = 3, force.exact = TRUE,
                     max.alpha = 4, max.digit = Inf,
                     phonetic = TRUE, encoding = c("primary", "alternate"),
                     phon.min.alpha = 5, min.enc = 3,
                     semantic = FALSE, syn = NULL) {
  # Preliminary Checks
  ###################################################################
  # Check method argument
  method <- match.arg(method)
  fields <- list(k1 = NULL, k2 = NULL)
  # Check excep argument
  if (!is.null(excep) && is.vector(excep, mode = "character") == FALSE) {
    stop('"excep" is not a character vector')
  }
  if (!is.null(excep)) {
    excep <- toupper(excep)
  } else {
    excep <- ""
  }
  # Check if kwic 1 is present
  if (is.null(kwic1)) {
    stop('"kwic1" is missing')
  }
  # Check if kwic1 is of class KWIC
  if (is(kwic1, "KWIC")) {
    fields[[1]] <- kwic1[[3]]
    kwic1 <- kwic1[[1]][!(kwic1[[1]]$KEYWORD %in% excep),c(1,3)]
  } else {
    stop('"kwic1" is not of class KWIC')
  }
  # Assign and check query and source kwic indexes according to methods
  if (method == "a") {
    kwic2 <- NULL
  }
  if (method == "b" | method == "c") {
    if (is.null(kwic2)) {
      stop('"kwic2" is missing')
    }
    if (is(kwic2, "KWIC")) {
      fields[[2]] <- kwic2[[3]]
      kwic2 <- kwic2[[1]][!(kwic2[[1]]$KEYWORD %in% excep),c(1,3)]
    } else {
      stop('"kwic2" is not of class KWIC')
    }
  }

  # Check encoding argument
  if (phonetic) {
    encoding <- match.arg(encoding)
  }
  # Check syn argument
  if (semantic) {
    if (is.null(syn)) {
      stop('"syn" is missing')
    }
    if (!is.list(syn)) {
      stop('"syn" is not a list')
    }
    if (is.element(FALSE, as.logical(lapply(syn,
                                           function(x) is.character(x))))) {
      warning('list "syn" had non character vectors; coerced to character')
      syn <-  as.logical(lapply(syn, function(x) as.character(x)))
    }
    nsyn1 <- length(syn)
    syn <- fix.syn(syn)
    nsyn2 <- length(syn)
    if (nsyn1 != nsyn2) {
      warning('synsets in list "syn" with common strings were merged')
    }
    rm(nsyn1, nsyn2)
  }
  # Create the output list
  out <- list(FuzzyDuplicates = NULL, PhoneticDuplicates = NULL,
              SemanticDuplicates = NULL)
  attr(out, "method") <- method
  attr(out, "fields") <- fields
  rm(fields)
  # Prepare the indexes
  ###################################################################
  kwicQ <- as.data.table(kwic1)
  rm(kwic1)
  kwicQ[, PRIM_ID := gsub("([[:space:]])\\1*", "", PRIM_ID)]
  kwicQ[, PRIM_ID := paste("[K1]", PRIM_ID, sep = "")]
  kwicQ[, IDKW := paste(PRIM_ID, KEYWORD, sep = ":")]
  if (method == "a" | method == "b") {
    kwicQ <- kwicQ[, list(PRIM_ID = paste0(setdiff(sort(unique(unlist(strsplit(get("PRIM_ID"),
                                                                               split = ", ")))), ""), collapse = ", "),
                          IDKW = paste0(setdiff(sort(unique(unlist(strsplit(get("IDKW"),
                                                                            split = ", ")))), ""), collapse = ", ")),
                   by = "KEYWORD"]
    # Id the chunks in kwic1
    M <- nrow(kwicQ)
    if (M > chunksize) {
      chunksize <- chunksize
    } else {
      chunksize <- M
    }
    kwicQ[, iter := rep(1:M, each = chunksize, length.out = M)]
    setcolorder(kwicQ, neworder = c("PRIM_ID", "KEYWORD", "IDKW", "iter"))
    setDF(kwicQ)
  }
  if (method == "a") {
    kwicS <- kwicQ
    #N <- M
  } else {
    kwicS <- as.data.table(kwic2)
    rm(kwic2)
    kwicS[, PRIM_ID := gsub("([[:space:]])\\1*", "", PRIM_ID)]
    kwicS[, PRIM_ID := paste("[K2]", PRIM_ID, sep = "")]
    kwicS[, IDKW := paste(PRIM_ID, KEYWORD, sep = ":")]
    if (method == "b") {
      kwicS <- kwicS[, list(PRIM_ID = paste0(setdiff(sort(unique(unlist(strsplit(get("PRIM_ID"),
                                                                                 split = ", ")))), ""), collapse = ", "),
                            IDKW = paste0(setdiff(sort(unique(unlist(strsplit(get("IDKW"),
                                                                              split = ", ")))), ""), collapse = ", ")),
                     by = "KEYWORD"]
    }
    N <- nrow(kwicS)
    #kwicS[, iter:= rep(1:M, each = chunksize, length.out = N)]
    setcolorder(kwicS, neworder = c("PRIM_ID", "KEYWORD", "IDKW"))
    setDF(kwicS)
    if (method == "c") {
      kwicQ <- as.data.table(rbind(kwicQ, kwicS))
      kwicQ <- kwicQ[, list(PRIM_ID = paste0(setdiff(sort(unique(unlist(strsplit(get("PRIM_ID"),
                                                                                 split= ", ")))), ""), collapse = ", "),
                            IDKW = paste0(setdiff(sort(unique(unlist(strsplit(get("IDKW"),
                                                                              split = ", ")))), ""), collapse = ", ")), by = "KEYWORD"]
      M <- nrow(kwicQ)
      if (M > chunksize) {
        chunksize <- chunksize
      } else {
        chunksize <- M
      }
      kwicQ[, iter := rep(1:M, each = chunksize, length.out = M)]
      setcolorder(kwicQ, neworder = c("PRIM_ID", "KEYWORD", "IDKW", "iter"))
      kwicS <- kwicQ
      setDF(kwicQ)
      setDF(kwicS)
    }
  }

  # Fuzzy Matching
  ###################################################################
  if (fuzzy) {
    # Coerce to integer
    if (max.alpha != Inf) {
      max.alpha <- as.integer(max.alpha)
    }
    if (max.digit != Inf) {
      max.digit <- as.integer(max.digit)
    }
    out[[1]] <- FuzzyDup(kwic1 = kwicQ, kwic2 = kwicS,
                         method = method, useBytes = useBytes,
                         max.dist = max.dist, force.exact = force.exact,
                         max.alpha = max.alpha, max.digit = max.digit)
  }
  # Phonetic Matching
  ###################################################################
  if (phonetic) {
    # Coerce to integer
    if (phon.min.alpha != Inf) {
      phon.min.alpha <- as.integer(phon.min.alpha)
    }
    if (min.enc != Inf) {
      min.enc <- as.integer(min.enc)
    }
    out[[2]] <- PhoneticDup(kwic1 = kwicQ, kwic2 = kwicS,
                            method = method, useBytes = useBytes,
                            encoding = encoding,
                            phon.min.alpha = phon.min.alpha, min.enc = min.enc)
  }
  # Semantic Matching
  ###################################################################
  if (semantic) {
    out[[3]] <- SemanticDup(kwic1 = kwicQ, kwic2 = kwicS,
                            method = method, syn = syn,
                            useBytes = useBytes)
  }
  class(out) <- "ProbDup"
  # Convert to null if no sets of a type are retrieved
  for (i in 1:3) {
    if (is.null(out[[i]])) {
      out[i] <- list(NULL)
    } else {
      if (dim(out[[i]])[1] == 0) {
        out[i] <- list(NULL)
      }
    }
  }
  return(out)
}

FuzzyDup <- function(kwic1, kwic2, max.dist, useBytes,
                     force.exact, max.alpha, max.digit, method) {
  kwic1 <- as.data.table(kwic1)
  kwic2 <- as.data.table(kwic2)
  M <- nrow(kwic1)
  N <- nrow(kwic2)
  # Identify strings to be exactly matched
  ind_exactQ <- logical(length = M)
  if (force.exact) {
    cond1 <- grepl("[[:digit:]]", kwic1$KEYWORD,
                   ignore.case = TRUE) == TRUE & stri_count_regex(kwic1$KEYWORD,
                                                                  "[[:digit:]]") <= max.digit
    cond2 <- grepl("[[:alpha:]]", kwic1$KEYWORD,
                   ignore.case = TRUE) == TRUE & stri_count_regex(kwic1$KEYWORD,
                                                                  "[[:alpha:]]") <= max.alpha
    ind_exactQ <- cond1 | cond2
    rm(cond1, cond2)
  }
  ind_exactS <- logical(length = N)
  if (force.exact) {
    cond1 <- grepl("[[:digit:]]", kwic2$KEYWORD,
                   ignore.case = TRUE) == TRUE & stri_count_regex(kwic2$KEYWORD,
                                                                  "[[:digit:]]") <= max.digit
    cond2 <- grepl("[[:alpha:]]", kwic2$KEYWORD,
                   ignore.case = TRUE) == TRUE & stri_count_regex(kwic2$KEYWORD,
                                                                  "[[:alpha:]]") <= max.alpha
    ind_exactS <- cond1 | cond2
    rm(cond1, cond2)
  }
  # Identify strings with mixed characters(alpha + digit)
  ind_mixedQ <- logical(length = M)
  cond1 <- grepl("[[:digit:]]", kwic1$KEYWORD,
                 ignore.case = TRUE) == TRUE & grepl("[[:alpha:]]",
                                                     kwic1$KEYWORD,
                                                     ignore.case = TRUE) == TRUE
  cond2 <- stri_count_regex(kwic1$KEYWORD, "[[:alpha:]]") > max.alpha
  ind_mixedQ <- cond1 & cond2
  rm(cond1, cond2)
  ind_mixedS <- logical(length = N)
  cond1 <- grepl("[[:digit:]]", kwic2$KEYWORD,
                 ignore.case = TRUE) == TRUE & grepl("[[:alpha:]]",
                                                     kwic2$KEYWORD,
                                                     ignore.case = TRUE) == TRUE
  cond2 <- stri_count_regex(kwic2$KEYWORD, "[[:alpha:]]") > max.alpha
  ind_mixedS <- cond1 & cond2
  rm(cond1, cond2)
  # Prepare mixed string vectors
  mixed_alphaQ <- ifelse(ind_mixedQ == TRUE, kwic1$KEYWORD, "")
  mixed_digitQ <- mixed_alphaQ
  mixed_alphaQ <- gsub(pattern = "[[:digit:]]", replacement = "",
                       x = mixed_alphaQ)
  mixed_digitQ <- gsub(pattern = "[[:alpha:]]", replacement = "",
                       x = mixed_digitQ)
  mixed_alphaS <- ifelse(ind_mixedS == TRUE, kwic2$KEYWORD, "")
  mixed_digitS <- mixed_alphaS
  mixed_alphaS <- gsub(pattern = "[[:digit:]]", replacement = "",
                       x = mixed_alphaS)
  mixed_digitS <- gsub(pattern = "[[:alpha:]]", replacement = "",
                       x = mixed_digitS)
  ind_mixed_digit_exact <- stri_count_regex(mixed_digitQ,
                                            "[[:digit:]]") <= max.digit
  # Create progress bar
  invisible(capture.output(pb <- txtProgressBar(min = 0, max = max(kwic1$iter),
                                                style = 3)))
  message("Fuzzy matching")
  for (i in unique(kwic1$iter)) {
    in_iter <- (kwic1$iter == i)
    # Compute the chunk stringdist matrices
    if (method == "b") {
      exact <- stringdistmatrix(a = kwic1$KEYWORD[in_iter & ind_exactQ],
                                b = kwic2$KEYWORD[ind_exactS],
                                method = "lv", useBytes = useBytes)
      exact[exact != 0] <- Inf
    }
    fuzzy <- stringdistmatrix(a = kwic1$KEYWORD[in_iter & !ind_exactQ],
                              b = kwic2$KEYWORD[!ind_exactS],
                              method = "lv", useBytes = useBytes)
    fuzzy[fuzzy > max.dist] <- Inf
    mixed <- NULL
    if (sum(in_iter & ind_mixedQ) > 0) {
      mixed_alpha_fuzzy <- stringdistmatrix(a = mixed_alphaQ[in_iter & ind_mixedQ],
                                            b = mixed_alphaS[ind_mixedS],
                                            method = "lv", useBytes = useBytes)
      mixed_alpha_fuzzy[mixed_alpha_fuzzy > max.dist] <- Inf
      mixed_digit_exact <- stringdistmatrix(a = mixed_digitQ[in_iter & ind_mixedQ & ind_mixed_digit_exact],
                                            b = mixed_digitS[ind_mixedS],
                                            method = "lv", useBytes = useBytes)
      mixed_digit_exact[mixed_digit_exact != 0] <- Inf
      mixed_digit_fuzzy <- stringdistmatrix(a = mixed_digitQ[in_iter & ind_mixedQ &  !ind_mixed_digit_exact],
                                            b = mixed_digitS[ind_mixedS],
                                            method = "lv", useBytes = useBytes)
      mixed_digit_fuzzy[mixed_digit_fuzzy > max.dist] <- Inf
#       if(sum(in_iter & ind_mixedQ) == 1) {
#         mixed_alpha_fuzzy <- t(mixed_alpha_fuzzy)
#       }
#       if(sum(in_iter & ind_mixedQ & ind_mixed_digit_exact) == 1) {
#         mixed_digit_exact <- t(mixed_digit_exact)
#       }
#       if(sum(in_iter & ind_mixedQ & !ind_mixed_digit_exact) == 1) {
#         mixed_digit_fuzzy <- t(mixed_digit_fuzzy)
#       }
      mixed_digit_comb <- matrix(numeric(), nrow = dim(mixed_alpha_fuzzy)[1],
                                 ncol = dim(mixed_alpha_fuzzy)[2])
      rownames(mixed_digit_comb) <- which(in_iter & ind_mixedQ)
      if (dim(mixed_digit_fuzzy)[1] != 0) {
        mixed_digit_comb[which(in_iter & ind_mixedQ &  !ind_mixed_digit_exact) %in% rownames(mixed_digit_comb),] <- mixed_digit_fuzzy
      }
      if (dim(mixed_digit_exact)[1] != 0) {
        mixed_digit_comb[which(in_iter & ind_mixedQ & ind_mixed_digit_exact) %in% rownames(mixed_digit_comb),] <- mixed_digit_exact
      }
      mixed <- mixed_alpha_fuzzy + mixed_digit_comb
      rm(mixed_alpha_fuzzy, mixed_digit_exact,
         mixed_digit_fuzzy, mixed_digit_comb)
    }
    # Checks
#     if (method == "b") {
#       if(sum(in_iter & ind_exactQ) == 1) {
#         exact <- t(exact)
#       }
#     }
#     if (sum(in_iter & !ind_exactQ) == 1) {
#       fuzzy <- t(fuzzy)
#     }
    # Fetch duplicates to a new column
    if (method == "b") {
      if (sum(in_iter & ind_exactQ ) > 0 & sum(ind_exactS) > 0) {
        kwic1[in_iter & ind_exactQ, FuzzydupIDKW := apply(exact, 1,
                                                         function(x) paste(as.character(unlist(kwic2[ind_exactS]$IDKW[x == 0])),
                                                                           collapse = ", "))]
      }
    }
    if (sum(in_iter & !ind_exactQ) > 0 & sum(!ind_exactS) > 0) {
      kwic1[in_iter & !ind_exactQ, FuzzydupIDKW := apply(fuzzy, 1,
                                                        function(x) paste(as.character(unlist(kwic2[!ind_exactS]$IDKW[x <= max.dist])),
                                                                          collapse = ", "))]
    }
    if (sum(in_iter & ind_mixedQ) > 0 & sum(ind_mixedS) > 0) {
      kwic1[in_iter & ind_mixedQ, FuzzydupIDKW2 := apply(mixed, 1,
                                                        function(x) paste(as.character(unlist(kwic2[ind_mixedS]$IDKW[x <= max.dist])),
                                                                          collapse = ", "))]
    }
  if (method == "b") {
    rm(exact)
    }
    rm(fuzzy, mixed)
    # update progress bar
    setTxtProgressBar(pb, i)
    cat("\rBlock", i, "/", max(kwic1$iter), "|")
  }
  close(pb)
  ind_exactQ2 <- stri_count_fixed(kwic1$PRIM_ID, ",") != 0 & ind_exactQ
  rm(ind_exactQ, ind_mixedQ, ind_exactS, ind_mixedS,
     pb, mixed_digitS, mixed_alphaS, mixed_alphaQ, mixed_digitQ,
     ind_mixed_digit_exact, in_iter)
  cols <- setdiff(colnames(kwic1), c("KEYWORD", "PRIM_ID", "iter"))
  for (j in cols) {
    set(kwic1,which(is.na(kwic1[[j]])),j,"")
  }
  if (method == "a" | method == "c") {
    kwic1[ind_exactQ2, FuzzydupIDKW := IDKW]
  }
  rm(ind_exactQ2)
  if ("FuzzydupIDKW2" %in% colnames(kwic1)) {
    kwic1[, FuzzydupIDKW := toString(unique(c(strsplit(FuzzydupIDKW,
                                                       split = ", ")[[1]],
                                              strsplit(FuzzydupIDKW2,
                                                       split = ", ")[[1]]))),
          by = IDKW]
    kwic1[, FuzzydupIDKW2 := NULL]
  }
  kwic1[, FuzzydupID := gsub(":\\S+\\b", "", FuzzydupIDKW)]
  kwic1 <- dupsets(kwic1, "F", method = method)
  return(kwic1)
}

PhoneticDup <- function(kwic1, kwic2, encoding, useBytes,
                        phon.min.alpha, min.enc, method) {
  kwic1 <- as.data.table(kwic1)
  kwic2 <- as.data.table(kwic2)
  M <- nrow(kwic1)
  N <- nrow(kwic2)
  # Fetch keyword phonetic encodings
  if (encoding == "primary") {
    DMQ <- DoubleMetaphone(kwic1$KEYWORD)[[1]]
    DMS <- DoubleMetaphone(kwic2$KEYWORD)[[1]]
  }
  if (encoding == "alternate") {
    DMQ <- DoubleMetaphone(kwic1$KEYWORD)[[2]]
    DMS <- DoubleMetaphone(kwic2$KEYWORD)[[2]]
  }
  # Prepare keyword phonetic encodings
  DMQ <- ifelse(stri_count_regex(kwic1$KEYWORD, "[[:alpha:]]") < phon.min.alpha,
                "", DMQ)
  DMS <- ifelse(stri_count_regex(kwic2$KEYWORD, "[[:alpha:]]") < phon.min.alpha,
                "", DMS)
  DMQ <- ifelse(nchar(DMQ) < min.enc, "", DMQ)
  DMS <- ifelse(nchar(DMS) < min.enc, "", DMS)
  # Identify strings with phonetic encodings
  ind_phon <- logical(length = M)
  ind_phon <- DMQ != ""
  ind_phonS <- logical(length = N)
  ind_phonS <- DMS != ""
  # Identify strings with phonetic encodings having digits
  ind_phon_digitQ <- grepl("[[:digit:]]", kwic1$KEYWORD,
                           ignore.case = TRUE) == TRUE & ind_phon
  ind_phon_digitS <- grepl("[[:digit:]]", kwic2$KEYWORD,
                           ignore.case = TRUE) == TRUE & ind_phonS
  # Prepared strings with phonetic encodings having digits
  phon_digitQ <- ifelse(ind_phon_digitQ == TRUE, kwic1$KEYWORD, "")
  phon_digitQ <- gsub(pattern = "[[:alpha:]]", replacement = "",
                      x = phon_digitQ)
  phon_digitS <- ifelse(ind_phon_digitS == TRUE, kwic2$KEYWORD, "")
  phon_digitS <- gsub(pattern = "[[:alpha:]]", replacement = "",
                      x = phon_digitS)
  # Create progress bar
  invisible(capture.output(pb <- txtProgressBar(min = 0, max = max(kwic1$iter),
                                                style = 3)))
  message("Phonetic matching")
  for (i in unique(kwic1$iter)) {
    in_iter <- (kwic1$iter == i)
    # Create distance matrix
    phon_dist <- stringdistmatrix(a = DMQ[in_iter & ind_phon],
                                  b = DMS[ind_phonS],
                                  method = "lv", useBytes = useBytes)
    phon_dist[phon_dist != 0] <- Inf
    # Checks
#     if(sum(in_iter & ind_phon) == 1) {
#       phon_dist <- t(phon_dist)
#     }
    if (sum(in_iter & ind_phon_digitQ) > 0) {
      phon_digit_exact <- stringdistmatrix(a = phon_digitQ[in_iter & ind_phon_digitQ],
                                           b = phon_digitS[ind_phon_digitS],
                                           method = "lv", useBytes = useBytes)
#       if(dim(phon_digit_exact)[2] == 1) {
#         phon_digit_exact <- t(phon_digit_exact)
#       }
      phon_digit_exact[phon_digit_exact != 0] <- Inf
      phon_dist_tr1 <-  which(which(in_iter & ind_phon) %in% which(in_iter & ind_phon_digitQ))
      phon_dist_tr2 <- which(which(ind_phonS) %in% which(ind_phon_digitS))
      phon_dist[phon_dist_tr1,phon_dist_tr2] <- phon_dist[phon_dist_tr1,phon_dist_tr2] + phon_digit_exact
      rm(phon_digit_exact)
    }
    # Fetch duplicates to a new column
    if (sum(in_iter & ind_phon) > 0 & sum(ind_phonS) > 0) {
      kwic1[in_iter & ind_phon, PhoneticdupIDKW := apply(phon_dist, 1,
                                                        function(x) paste(as.character(unlist(kwic2[ind_phonS]$IDKW[x == 0])),
                                                                          collapse = ", "))]
    }
    rm(phon_dist,in_iter)
    # update progress bar
    setTxtProgressBar(pb, i)
    cat("\rBlock", i, "/", max(kwic1$iter), "|")
  }
  close(pb)
  rm(ind_phon, pb, ind_phonS, ind_phon_digitQ, phon_digitQ, phon_digitS)
  for (j in c("IDKW", "PhoneticdupIDKW")) {
    set(kwic1,which(is.na(kwic1[[j]])),j,"")
  }
  kwic1[, PhoneticdupID := gsub(":\\S+\\b", "", PhoneticdupIDKW)]
  kwic1 <- dupsets(kwic1, "P", method = method)
  return(kwic1)
}

SemanticDup <- function(kwic1, kwic2, syn, useBytes, method) {
  kwic1 <- as.data.table(kwic1)
  kwic2 <- as.data.table(kwic2)
  M <- nrow(kwic1)
  N <- nrow(kwic2)
  # Identify keyword strings associated with synsets
  SMQ <- as.character(with(stack(syn), ind[match(kwic1$KEYWORD, values)]))
  SMS <- as.character(with(stack(syn), ind[match(kwic2$KEYWORD, values)]))
  SMQ[is.na(SMQ)]   <- ""
  SMS[is.na(SMS)]   <- ""
  # Identify strings for semantic matching
  ind_sem <- SMQ != ""
  ind_semS <- SMS != ""
  # Create progress bar
  invisible(capture.output(pb <- txtProgressBar(min = 0, max = max(kwic1$iter),
                                                style = 3)))
  message("Semantic matching")
  for (i in unique(kwic1$iter)) {
    in_iter <- (kwic1$iter == i)
    # Create distance matrix
    sem_dist <- stringdistmatrix(a = SMQ[in_iter & ind_sem],
                                 b = SMS[ind_semS],
                                 method = "lv", useBytes = useBytes)
    sem_dist[sem_dist != 0] <- Inf
    # Checks
#     if (sum(in_iter & ind_sem) == 1) {
#       sem_dist <- t(sem_dist)
#     }
    # Fetch duplicates to a new column
    if (sum(in_iter & ind_sem) > 0 & sum(ind_semS) > 0) {
      kwic1[in_iter & ind_sem, SemanticdupIDKW := apply(sem_dist, 1,
                                                        function(x) paste(as.character(unlist(kwic2[ind_semS]$IDKW[x == 0])),
                                                                          collapse = ", "))]
    }
    rm(sem_dist,in_iter)
    # update progress bar
    setTxtProgressBar(pb, i)
    cat("\rBlock", i, "/", max(kwic1$iter), "|")
  }
  close(pb)
  rm(ind_sem, ind_semS)
  for (j in c("IDKW", "SemanticdupIDKW")) {
    set(kwic1,which(is.na(kwic1[[j]])),j,"")
  }
  kwic1[, SemanticdupID := gsub(":\\S+\\b", "", SemanticdupIDKW)]
  kwic1 <- dupsets(kwic1, "S", method = method)
  # Remove synsets with single/unique members
  kwic1 <- kwic1[!unlist(lapply(strsplit( gsub("*?\\[\\S+:", "", kwic1$IDKW),
                                          split = ", "), function(x) length(unique(x)))) == 1,]
  return(kwic1)
}

dupsets <- function(kwicout, type, method) {
  if (dim(kwicout)[2] == 6) {
    kwicout <- as.data.table(subset(kwicout, get(names(kwicout)[6]) != ""))
    setkey(kwicout, PRIM_ID)
    kwicout[, c("iter", "KEYWORD") := NULL]
    if (method == "b") {
      kwicout[,  3 := paste(get(names(kwicout)[2]),
                            get(names(kwicout)[3]), sep = ", "), with = FALSE]
      kwicout[,  4 := paste(get(names(kwicout)[1]),
                            get(names(kwicout)[4]), sep = ", "), with = FALSE]
    }
    kwicout[, IDKW := NULL]
    if (method != "b") {
      kwicout[, Ndup := stri_count_fixed(get(colnames(kwicout)[3]), ",")]
      kwicout <- subset(kwicout, Ndup != 0)
      kwicout[, Ndup := NULL]
    }
    # Merge by PRIM_ID, then by ID, then add TYPE
    kwicout <- kwicout[, list(ID = paste0(setdiff(sort(unique(unlist(strsplit(get(names(kwicout)[3]), split = ", ")))), ""),
                                          collapse = ", "),
                              IDKW = paste0(setdiff(sort(unique(unlist(strsplit(get(names(kwicout)[2]), split = ", ")))), ""),
                                            collapse = ", ")),
                       by = "PRIM_ID"][, list(IDKW = paste0(sort(unique(unlist(strsplit(IDKW, split = ", ")))),
                                                            collapse = ", "),
                                                     TYPE = type), by = "ID"]
    setkey(kwicout, NULL)
    kwicout <- unique(kwicout)
    # Add SET_NO
    setkey(kwicout, "ID")
    kwicout[, SET_NO := as.factor(ID)]
    kwicout[, SET_NO := as.numeric(SET_NO)]
    setkey(kwicout, "SET_NO")
    # Add count
    kwicout[, COUNT := stri_count_fixed(ID, ",") + 1]
    kwicout <- subset(kwicout, COUNT != 1)
    # Finalise output
    setcolorder(kwicout, c("SET_NO", "TYPE", "ID", "IDKW", "COUNT"))
    setkey(kwicout, "SET_NO")
    setDF(kwicout)
    return(kwicout)
  }
}


fix.syn <- function(syn) {
  names(syn) <- NULL
  names(syn) <- paste0("SM", seq_along(syn))
  syn <- lapply(syn, sort)
  syn <- subset(syn, duplicated(syn) == FALSE)
  if (is.element(TRUE, sapply(syn, function(x) length(x) == 1))) {
    syn <- syn[!sapply(syn, function(x) length(x) == 1)]
    warning("synsets encountered in list 'syn' with length 1")
  }
  syncomb <- do.call("rbind",lapply(syn, embed, 2))
  gg <- graph.edgelist(syncomb, directed = F)
  x <- split(V(gg)$name, clusters(gg)$membership)
  x <- lapply(x, sort)
  x <- subset(x, duplicated(x) == FALSE)
  names(x) <- paste0("SM", seq_along(x))
  return(x)
}
