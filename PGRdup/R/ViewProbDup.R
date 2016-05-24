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

#' Visualize the probable duplicate sets retrieved in a \code{ProbDup} object
#' 
#' \code{ViewProbDup} plots summary visualizations of accessions within the
#' probable duplicate sets retrieved in a \code{ProbDup} object according to a 
#' grouping factor field(column) in the original database(s).
#' 
#' @param pdup An object of class \code{ProbDup}.
#' @param db1 A data frame of the PGR passport database.
#' @param db2 A data frame of the PGR passport database. Required when 
#'   \code{pdup} was created using more than one KWIC Index.
#' @param factor.db1 The \code{db1} column to be considered for grouping the 
#'   accessions. Should be of class character or factor.
#' @param factor.db2 The \code{db2} column to be considered for grouping the 
#'   accessions. Should be of class character or factor. retrieved.
#' @param max.count The maximum count of probable duplicate sets whose 
#'   information is to be plotted (see \strong{Note}).
#' @param select A character vector of factor names in \code{factor.db1} and/or 
#'   \code{factor.db2} to be considered for grouping accessions (see 
#'   \strong{Note}).
#' @param order The order of the type of sets retrieved in the plot. The default
#'   is \code{"type"} (see \strong{Details}).
#' @param main The title of the plot.
#'   
#' @return A list containing the following objects: \tabular{ll}{ 
#'   \code{Summary1} \tab The summary \code{data.frame} of number of accessions 
#'   per factor level. \cr \code{Summary2} \tab The summary \code{data.frame} of
#'   number of accessions and sets per each type of sets classified according to
#'   factor levels. \cr \code{SummaryGrob} \tab A grid graphical object (Grob) 
#'   of the summary visualization plot. \cr }
#'   
#' @note When any primary ID/key records in the fuzzy, phonetic or semantic 
#'   duplicate sets are found to be missing from the original databases 
#'   \code{db1} and \code{db2}, then they are ignored and only the matching 
#'   records are considered for visualization.
#'   
#'   This may be due to data standardization of the primary ID/key field using 
#'   the function \code{\link[PGRdup]{DataClean}} before creation of the KWIC 
#'   index and subsequent identification of probable duplicate sets. In such a 
#'   case, it is recommended to use an identical data standardization operation 
#'   on the databases \code{db1} and \code{db2} before running this function. 
#'   For summary and visualization of the set information in the object of class
#'   \code{ProbDup} by \code{ViewProbDup}, the disjoint of the retrieved sets 
#'   are made use of, as they are more meaningful than the raw sets retrieved. 
#'   So it is recommended that the disjoint of sets obtained using the 
#'   \code{DisProbDup} be used as the input \code{pdup}.
#'   
#'   All the accession records in sets with count > \code{max.count} will be 
#'   considered as being unique.
#'   
#'   The factor levels in the \code{factor.db1} and/or \code{factor.db2} columns
#'   corresponding to those mentioned in \code{select} argument alone will be 
#'   considered for visualization. All other factor levels will be grouped 
#'   together to a single level named "Others".
#'   
#'   The argument \code{order} can be used to specify the order in which the 
#'   type of sets retrieved are to be plotted in the visualization. The default 
#'   \code{"type"} will order according to the kind of sets, \code{"sets"} will 
#'   order according to the number of sets in each kind and \code{"acc"} will 
#'   order according to the number of accessions in each kind.
#'   
#'   The individual plots are made using \code{\link[ggplot2]{ggplot}} and then 
#'   grouped together using \code{\link[gridExtra]{gridExtra-package}}.
#'   
#' @examples
#' \dontrun{
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
#' GN1$SourceCountry <- toupper(GN1$SourceCountry)
#' GN2$SourceCountry <- toupper(GN2$SourceCountry)
#' 
#' GN1$SourceCountry <- gsub("UNITED STATES OF AMERICA", "USA", GN1$SourceCountry)
#' GN2$SourceCountry <- gsub("UNITED STATES OF AMERICA", "USA", GN2$SourceCountry)
#' 
#' # Specify as a vector the database fields to be used
#' GN1fields <- c("NationalID", "CollNo", "OtherID1", "OtherID2")
#' GN2fields <- c("DonorID", "CollNo", "OtherID1", "OtherID2")
#' 
#' # Clean the data
#' GN1[GN1fields] <- lapply(GN1[GN1fields], function(x) DataClean(x))
#' GN2[GN2fields] <- lapply(GN2[GN2fields], function(x) DataClean(x))
#' y1 <- list(c("Gujarat", "Dwarf"), c("Castle", "Cary"), c("Small", "Japan"),
#'            c("Big", "Japan"), c("Mani", "Blanco"), c("Uganda", "Erect"),
#'            c("Mota", "Company"))
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
#'           "DARK", "E", "EARLY", "EC", "ERECT", "EXOTIC", "FLESH", "GROUNDNUT",
#'           "GUTHUKAI", "IMPROVED", "K", "KUTHUKADAL", "KUTHUKAI", "LARGE",
#'           "LIGHT", "LOCAL", "OF", "OVERO", "P", "PEANUT", "PURPLE", "R",
#'           "RED", "RUNNER", "S1", "SAM", "SMALL", "SPANISH", "TAN", "TYPE",
#'           "U", "VALENCIA", "VIRGINIA", "WHITE")
#' 
#' # Specify the synsets as a list
#' syn <- list(c("CHANDRA", "AH114"), c("TG1", "VIKRAM"))
#' 
#' GNdupc <- ProbDup(kwic1 = GN1KWIC, kwic2 = GN2KWIC, method = "c",
#'                   excep = exep, fuzzy = TRUE, phonetic = TRUE,
#'                   encoding = "primary", semantic = TRUE, syn = syn)
#' 
#' GNdupcView <- ViewProbDup(GNdupc, GN1, GN2, "SourceCountry", "SourceCountry",
#'                          max.count = 30, select = c("INDIA", "USA"), order = "type",
#'                          main = "Groundnut Probable Duplicates")
#' 
#' }       
#' @seealso \code{\link[PGRdup]{ProbDup}}, \code{\link[PGRdup]{DisProbDup}}, 
#'   \code{\link[PGRdup]{DataClean}}, \code{\link[ggplot2]{ggplot}}, 
#'   \code{\link[gridExtra]{gridExtra-package}}
#'   
#' @import data.table
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ggplotGrob
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 coord_polar
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 labs
#' @importFrom grid textGrob
#' @importFrom grid unit
#' @importFrom gridExtra arrangeGrob
#' @importFrom methods is
#' @importFrom utils head
#' @export


ViewProbDup <- function(pdup, db1, db2 = NULL,
                             factor.db1, factor.db2 = NULL,
                             max.count = 30, select, order = "type",
                             main = NULL) {
  # Preliminary Checks
  if (!is(pdup, "ProbDup")) {
    stop("'pdup' is not of class ProbDup")
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
  if (!is.null(factor.db1) && !is.vector(factor.db1, mode = "character")) {
    stop("'factor.db1' is not a character vector")
  }
  if (is.element(FALSE, factor.db1 %in% colnames(db1))) {
    # Check if factor fields are present in db1 and stop if not
    stop("The 'factor.db1' field is missing in 'db1'.")
  }
  if (!(is.character(db1[, factor.db1]) || is.factor(db1[, factor.db1]))) {
    # check if factor column vector is of type character/vector
    stop("The 'factor.db1' column is not of type character or factor")
  }
  if (FALSE %in% (select %in% db1[, factor.db1])) {
    Pt <- data.frame(select, check = select %in% db1[, factor.db1], stringsAsFactors = F)
    stop(paste("The following selected factor(s) is/are missing from 'factor.db1' column",
               paste(Pt[Pt$check == FALSE, 1], collapse = ", "), sep = "\n"))
  }
  fields[[1]] <- union(fields[[1]], factor.db1)
  #setDT(db1)
  db1 <- as.data.table(db1)
  if (is.factor(db1[, get(factor.db1)])) {
    db1[, (factor.db1) := as.character(get(factor.db1))]
  }
  # Convert whitespace to NA
  db1[grepl("^[[:space:]]*$", get(factor.db1)), (factor.db1) := NA]
  if (!identical(setdiff(colnames(db1), fields[[1]]), character(0))) {
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
      stop('One or more kwic2 fields are missing in "db2"')
    }
    if (!is.null(factor.db2) && !is.vector(factor.db2, mode = "character")) {
      stop('"factor.db2" is not a character vector')
    }
    if (is.element(FALSE, factor.db2 %in% colnames(db2)) == TRUE) {
      # Check if factor fields are present in db2 and stop if not
      stop('The "factor.db2" field is missing in "db2"')
    }
    if (!(is.character(db2[, factor.db2]) || is.factor(db2[, factor.db2]))) {
      # check if factor column vector is of type character/vector
      stop("The 'factor.db2' column is not of type character or factor")
    }
    if (FALSE %in% (select %in% db2[, factor.db2])) {
      Pt <- data.frame(select, check = select %in% db2[, factor.db2],
                       stringsAsFactors = F)
      stop(paste("The following selected factor(s) is/are missing from 'factor.db2' column",
                 paste(Pt[Pt$check == FALSE, 1], collapse = ", "), sep = "\n"))
    }

    fields[[2]] <- union(fields[[2]], factor.db2)
    #setDT(db2)
    db2 <- as.data.table(db2)
    if (is.factor(db2[, get(factor.db2)])) {
      db2[, (factor.db2) := as.character(get(factor.db2))]
    }
    # Convert whitespace to NA
    db2[grepl("^[[:space:]]*$", get(factor.db2)), (factor.db2) := NA]
    if (!identical(setdiff(colnames(db2), fields[[1]]), character(0))) {
      db2[, setdiff(colnames(db2), fields[[2]]) := NULL]
    }
    db2[, K2_PRIM_ID := get(fields[[2]][1])]
    setcolorder(db2, neworder = union("K2_PRIM_ID",
                                      setdiff(colnames(db2), "K2_PRIM_ID")))
    db2[, K2_PRIM_ID := as.character(K2_PRIM_ID)]
    setkey(db2, "K2_PRIM_ID")
  }
  order <- match.arg(order, c("type", "sets", "acc"), several.ok = FALSE)
  # Get disjoint sets
  if (!is.null(pdup[["DisjointDupicates"]])) {
    tryCatch(DisProbDup(pdup),
             warning = function(e) {
               message(gsub("returned", "displayed", e$message))
             })
    pdup <- suppressWarnings(DisProbDup(pdup))
  } else {
    tstr <- data.frame(names = c("F", "P", "S"),
                       `is.null` = c(is.null(pdup[[1]]), is.null(pdup[[2]]),
                                     is.null(pdup[[3]])),
                       stringsAsFactors = F)
    if (dim(tstr[tstr$is.null == FALSE, ])[1] == 1) {
      tryCatch(DisProbDup(pdup),
      warning = function(e) {
        message(gsub("returned", "displayed", e$message))
      })
      pdup <- suppressWarnings(DisProbDup(pdup))
    } else {
      pdup <- DisProbDup(pdup, combine = tstr[tstr$is.null == FALSE,]$names)
    }
  }

  # Get Review data.frame
  tryCatch(ReviewProbDup(pdup, db1, db2,
                         extra.db1 = factor.db1, extra.db2 = factor.db2,
                         max.count, insert.blanks = FALSE),
           warning = function(e) {
             message(gsub("merged", "used", e$message))
           })
  sets <- suppressWarnings(ReviewProbDup(pdup, db1, db2,
                                         extra.db1 = factor.db1,
                                         extra.db2 = factor.db2,
                                         max.count, insert.blanks = FALSE))
  # Prepare the sets
  if (is.null(factor.db2)) {
    cols <- c("SET_NO", "PRIM_ID",
              paste("K1X", factor.db1, sep = "_"))
  } else {
    cols <- c("SET_NO", "PRIM_ID",
              paste("K1X", factor.db1, sep = "_"),
              paste("K2X", factor.db2, sep = "_"))
  }
  sets <- as.data.table(sets[, cols])
  #setDT(sets)
  if (is.null(factor.db2)) {
    sets[, FACTOR := get(cols[3])]
  } else {
    sets[, FACTOR := ifelse(is.na(get(cols[3])), get(cols[4]), get(cols[3]))]
  }
  sets[!FACTOR %in% select, FACTOR := "Others"]
  sets[, FACTOR := as.factor(FACTOR)]
  sets <- sets[order(SET_NO), as.list(table(FACTOR)),
               by = SET_NO] # Get summary table
  factorlevels <- setdiff(colnames(sets), "SET_NO")
  Nfactorlevels <- paste0("No.Acc. (", factorlevels, ")")
  setnames(sets, old = factorlevels, new = Nfactorlevels)
  sets[, (factorlevels) := data.frame(matrix(ifelse(sets[,Nfactorlevels, with = FALSE] != 0, 1, 0),
                                             ncol = length(Nfactorlevels)))]
  sets[, WITHIN := ifelse(rowSums(sets[, factorlevels, with = FALSE]) == 1,
                          1, 0)]
  sets[, `NO.ACC` := rowSums(sets[, Nfactorlevels, with = FALSE])]
  # Prepare orphans
  dups <- ParseProbDup(pdup, insert.blanks = F, max.count)$PRIM_ID
  orps <- db1[!K1_PRIM_ID %in% dups, .SD, .SDcols = c("K1_PRIM_ID",
                                                      factor.db1)]
  setnames(orps, colnames(orps), c("PRIM_ID", "FACTOR"))
  if (!is.null(factor.db2)) {
    orp2 <- db2[!K2_PRIM_ID %in% dups, .SD, .SDcols = c("K2_PRIM_ID",
                                                        factor.db2)]
    setnames(orp2, colnames(orp2), c("PRIM_ID", "FACTOR"))
    orps <- rbind(orps, orp2)
    rm(orp2)
  }
  ndupsets <- sum(unlist(sapply(pdup, function(x) dim(x)[1]))) # Add set numbers
  orps[, SET_NO := 1:nrow(orps) + ndupsets]
  orps[!FACTOR %in% select, FACTOR := "Others"]
  orps[, FACTOR := as.factor(FACTOR)]
  orps <- orps[order(SET_NO), as.list(table(FACTOR)),
               by = SET_NO] # Get summary table
  factorlevels <- setdiff(colnames(orps), "SET_NO")
  Nfactorlevels <- paste0("No.Acc. (", factorlevels, ")")
  orps[, WITHIN := 0][, `NO.ACC` := 1]
  orps[, (Nfactorlevels) := orps[, factorlevels, with = FALSE]]
  # Final summary = sets + orphans
  bby <- append(factorlevels, "WITHIN")
  setcolorder(sets, colnames(orps))
  sets <- rbind(sets, orps)
  setkeyv(sets, append(factorlevels, "WITHIN"))
  summ1 <- sets[, lapply(.SD, sum), .SDcols = Nfactorlevels, by = bby]
  summ2 <- sets[, list(NO.ACC = sum(NO.ACC),NO.SETS = .N), by = bby]
  summ <- merge(summ1, summ2)
  summ <- summ[, lapply(.SD, as.numeric)]
  summ[, KIND := as.factor(rownames(summ))]
  rm(summ1, summ2)
  # Specify the shape for mat plot
  summ[, singles := ifelse(rowSums(summ[, factorlevels, with = FALSE]) == 1,
                           1, 0)]
  summ[, shape :=  singles + WITHIN]
  # Reorder summary table
  if (order == "type") {
    tord <- as.numeric(append(rev(as.character(summ[summ$shape == 0]$KIND)),
                              append(as.character(summ[summ$shape == 2]$KIND),
                                     as.character(summ[summ$shape == 1]$KIND))))
    summ[, KIND := factor(KIND, levels = summ[tord, KIND])]
  }
  if (order == "acc") {
    summ[, KIND := factor(KIND, levels = summ[order(NO.ACC), KIND])]
  }
  if (order == "sets") {
    summ[, KIND := factor(KIND,
                          levels = summ[order(NO.SETS, na.last = FALSE), KIND])]
  }
  # Set blanks in No.sets plot
  summ[which(rowSums(summ[, factorlevels, with = F]) == 1),
       NO.SETS := ifelse(WITHIN == 0, yes = NA, NO.SETS)]
  # Create the data.frame to specify shading
  shading <- data.frame(min = seq(from = 0.5,
                                  to = max(as.numeric(summ$KIND)),
                                  by = 1),
                        max = seq(from = 1.5,
                                  to = max(as.numeric(summ$KIND)) + 0.5,
                                  by = 1),
                        col = 0)
  shading[which(seq(1:nrow(shading)) %% 2 == 1), ]$col <- 1
  shading$col <- ifelse(shading$col, "gray83", "gray93")
  #Stacked bar
  sbar <- melt(summ, id.vars = "KIND", measure.vars = Nfactorlevels)
  sbarg <- ggplot() +
    theme(panel.background = element_rect(fill = "transparent"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), #<>#
          axis.ticks.y = element_line(colour = "black"),
          axis.line = element_line(colour = NA),
          axis.line.y = element_line(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major.x = element_line(colour = "gray83"),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0,0.5,0.1,0.5), "cm")) +
    geom_bar(sbar, mapping = aes(x = KIND, y = value, fill = variable),
             stat = "identity") +
    geom_rect(data = shading,
              aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, alpha = 0.1),
              fill = shading$col, colour = "white") +
    geom_bar(sbar, mapping = aes(x = KIND, y = value, fill = variable),
             stat = "identity") +
    ylab("No. of accessions") +
    #coord_flip() +
    guides(fill = FALSE, colour = FALSE, alpha = FALSE)
  #dot plot
  dotp <- melt(summ, id.vars = "KIND", measure.vars = "NO.SETS")
  dotpg <- ggplot() +
    theme(panel.background = element_rect(fill = "transparent"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), #<>#
          axis.line = element_line(colour = NA),
          axis.line.y = element_line(colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major.x = element_line(colour = "gray83"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "gray33",
                                            linetype = "dotted"),
          plot.margin = unit(c(1,0.5,0.3,0.5), "cm")) +
    geom_point(dotp, mapping = aes(x = KIND, y = value), colour = "gray23",
               size = 3, pch = 18, na.rm = TRUE) +
    geom_rect(data = shading,
              aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, alpha = 0.1),
              fill = shading$col, colour = "white") +
    geom_point(dotp, mapping = aes(x = KIND, y = value), colour = "gray23",
               size = 3, pch = 18, fill = "gray23", na.rm = TRUE) +
    ylab("No. of sets") +
    #coord_flip() +
    #scale_y_reverse() +
    guides(fill = FALSE, colour = FALSE, alpha = FALSE)
  # matrix
  mat <- melt(summ, id.vars = c("KIND", "WITHIN", "shape", "singles"),
              measure.vars = factorlevels)
  mat[, shape := ifelse(shape == 0, "Between",
                        ifelse(shape == 1, "Unique", "Within"))]
  mat[, shape2 := ifelse(value == 1, shape, NA)]
  mat[, KIND2 := as.numeric(ifelse(is.na(shape2), NA, KIND))]
  mat[is.na(mat$KIND2), KIND2 := seq(from = max(mat$KIND2, na.rm = T) + 1,
                                     length.out = length(mat[is.na(KIND2),]$KIND2))]
  matg <- ggplot() +
    theme(panel.background = element_rect(fill = "transparent"),
          axis.title.x =  element_text(colour = "black",
                                       margin = margin(20,0,20,0)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), #<>#
          axis.ticks.y = element_blank(),
          axis.line = element_line(colour = NA),
          axis.line.y = element_line(colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_blank(), #<>#
          #legend.position = "none", #<>#
          plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
          #legend.title = element_blank(),
          panel.grid = element_blank()) +
    geom_tile(data = mat, mapping = aes(y = variable, x = as.factor(KIND)),
              colour = "white", fill = NA) +
    geom_rect(data = shading,
              aes(xmin = min, xmax = max,
                  ymin = 0.5, ymax = length(factorlevels) + 0.5, alpha = 0.1),
              fill = shading$col, colour = "white") +
    geom_tile(data = mat, mapping = aes(y = variable, x = as.factor(KIND)),
              colour = "white", fill = NA) +
    geom_point(data = mat, aes(y = variable, x = KIND, shape = mat$shape),
               colour = "gray73", alpha = 0.6, size = 5, na.rm = TRUE) +
    geom_point(data = mat, aes(y = variable, x = KIND, shape = mat$shape2),
               colour = "gray23",  size = 5,  na.rm = TRUE) +
    geom_line(data = mat, aes(y = variable, x = KIND, shape = mat$shape2,
                              group = KIND2),
              colour = "gray23",  size = 1,  na.rm = TRUE) +
    #coord_flip() +
    guides(shape = guide_legend(title = "Accession\nType"),
           fill = FALSE, colour = FALSE, alpha = FALSE) +
    xlab("Type of sets retrieved")
  # Extract legend1
  grobs <- ggplotGrob(matg)$grobs
  legend1 <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  rm(grobs)
  # Accession summary
  nacc <- colSums(summ[, Nfactorlevels, with = FALSE])
  nacc <- data.frame(nacc)
  nacc$factor <- rownames(nacc)
  nacc$factor <- factorlevels
  naccg <- ggplot(nacc, aes(x = factor, y = nacc, fill = factor)) +
    theme(panel.background = element_rect(fill = "transparent"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_line(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black",
                                      margin = margin(20,0,20,0)),
          panel.grid.major.y = element_line(colour = "gray93"),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0,0,0.5,0), "cm")) +
    geom_bar(stat = "identity", position = "identity") +
    ylab("No. of accessions") +
    coord_flip() +
    guides(fill  = guide_legend(title = "Accession\nGroup"),
           colour = FALSE, alpha = FALSE) +
    scale_y_reverse()
  # Extract legend1
  grobs <- ggplotGrob(naccg)$grobs
  legend2 <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  rm(grobs)
  # Summary donut
  smry <- data.frame(Total = c("Sets", "Accessions", "  Redundant", "  Unique"),
                     Count = c(sum(dotp[!is.na(dotp$value), ]$value),
                               sum(nacc$nacc),
                               sum(sbar[sbar$KIND %in% dotp[!is.na(dotp$value), ]$KIND, ]$value),
                               sum(sbar[!sbar$KIND %in% dotp[!is.na(dotp$value), ]$KIND, ]$value)))
  smry1 <- smry[3:4,]
  smry1$fraction <- smry1$Count / sum(smry1$Count)
  smry1$ymax <- cumsum(smry1$fraction)
  smry1$ymin <- c(0, head(smry1$ymax, n = -1))
  smry1$label <- paste(smry1$Total, "\n", round(smry1$fraction * 100, 2),
                       " %", sep = "")
  smryg <- ggplot(data = smry1, aes(fill = Total, ymax = ymax, ymin = ymin,
                                    xmax = 4, xmin = 3)) +
    geom_rect(show.legend = FALSE, fill = c("grey93", "grey73")) +
    coord_polar(theta = "y") +
    xlim(c(0, 4)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    geom_text(aes(x = 3.5, y = ( (ymin + ymax) / 2), label = label)) +
    xlab("") +
    ylab("") +
    labs(title = paste("\n\nTotal no. of sets = ", smry[1,]$Count, "\n",
                       "Total no. of accessions  = ", smry[2,]$Count,
                       sep = ""))
  # Get and modify grobs
  naccg1 <- ggplotGrob(naccg + theme(legend.position = "none"))
  dotpg1 <- ggplotGrob(dotpg)
  matg1 <- ggplotGrob(matg + theme(legend.position = "none"))
  sbarg1 <- ggplotGrob(sbarg)
  smryg1 <- ggplotGrob(smryg)
  smryg1$layout$clip[smryg1$layout$name == "panel"] <- "off"
  tempgrob <- cbind(naccg1, matg1, size = "last")
  naccg1$heights <- tempgrob$heights
  matg1$heights <- tempgrob$heights
  dotpg1$widths <- matg1$widths
  sbarg1$widths <- matg1$widths
  legend <-  arrangeGrob(legend1, legend2, ncol = 2)
  SummaryGrob <- arrangeGrob(smryg1, dotpg1,
                             legend, sbarg1,
                             naccg1, matg1,
                             ncol = 2, widths = c(0.5, 1.5),
                             heights = c(0.8,1,0.6),
                             top = textGrob(main, just = "top"))
  Summary1 <- data.frame("Factor" = nacc$factor,
                         "No. of accessions" = nacc$nacc,
                        check.names = FALSE, stringsAsFactors = FALSE)
  Summary2 <- as.data.table(summ)
  Summary2[, c("singles","shape") := NULL]
  setcolorder(Summary2, union("KIND", setdiff(colnames(Summary2), "KIND")))
  setnames(Summary2, old = "KIND", "Type of sets retrieved")
  setnames(Summary2, old = "NO.ACC", "Total no. of accessions")
  setnames(Summary2, old = "NO.SETS", "Total no. of sets")
  setDF(Summary2)
  out <- list(Summary1 = Summary1, Summary2 = Summary2,
              SummaryGrob = SummaryGrob)
  rm(Summary1, Summary2, SummaryGrob)
  return(out)
}
