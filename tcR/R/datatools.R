########## Support functions for managing the data ##########


#' Fix alleles / genes by removing allele information / unnecessary colons. 
#' 
#' @aliases fix.alleles fix.genes
#' 
#' @description
#' Fix alleles / genes by removing allele information / unnecessary colons. 
#' 
#' @param .data tcR data frame.
fix.alleles <- function (.data) {
  if (has.class(.data, "list")) {
    return(lapply(.data, fix.alleles))
  }
  
  .data$V.gene <- gsub("[*][[:digit:]]*", "", .data$V.gene)
  .data$D.gene <- gsub("[*][[:digit:]]*", "", .data$D.gene)
  .data$J.gene <- gsub("[*][[:digit:]]*", "", .data$J.gene)
  .data
}

fix.genes <- function (.data) {
  if (has.class(.data, "list")) {
    return(lapply(.data, fix.genes))
  }
  
  .fix <- function (.col) {
    # it's not a mistake
    .col <- gsub(", ", ",", .col, fixed = T, useBytes = T)
    .col <- gsub(",", ", ", .col, fixed = T, useBytes = T)
    .col
  }
  
  .data$V.gene <- .fix(.data$V.gene)
  .data$D.gene <- .fix(.data$D.gene)
  .data$J.gene <- .fix(.data$J.gene)
  .data
}


#' Print the given message if second parameter is a TRUE.
#' 
#' @param .message Character vector standing for a message.
#' @param .verbose If T then print the given mesasge.
#' @return Nothing.
.verbose.msg <- function (.message, .verbose = T) {
  if (.verbose) cat(.message)
}


#' Choose the right column.
#' 
#' @param .x Character vector with column IDs.
#' @param .verbose If T then print the error mesasge.
#' @return Character.
.column.choice <- function (x, .verbose = T) {
  x <- switch(x[1],
              read.count = "Read.count",
              umi.count = "Umi.count",
              read.prop = "Read.proportion",
              umi.prop = "Umi.proportion",
              { .verbose.msg("You have specified an invalid column identifier. Choosed column: Read.count\n", .verbose); "Read.count" })
  x
}


#' Fix names in lists.
#' 
#' @param .datalist List with data frames.
#' @return List with fixed names.
.fix.listnames <- function (.datalist) {
  if (is.null(names(.datalist))) { names(.datalist) <- paste0("Sample.", 1:length(.datalist)) }
  else {
    for (i in 1:length(.datalist)) {
      if (names(.datalist)[i] == "" || is.na(names(.datalist)[i])) {
        names(.datalist)[i] <- paste0("Sample.", i)
      }
    }
  }
  .datalist
}


#' Get all unique clonotypes.
#' 
#' @description
#' Get all unique clonotypes with merged counts. Unique clonotypes are those with
#' either equal CDR3 sequence or with equal CDR3 sequence and equal gene segments.
#' Counts of equal clonotypes will be summed up.
#' 
#' @param .data Either tcR data frame or a list with data frames.
#' @param .gene.col Either name of the column with gene segments used to compare clonotypes
#' or NA if you don't need comparing using gene segments.
#' @param .count.col Name of the column with counts for each clonotype.
#' @param .prop.col Name of the column with proportions for each clonotype.
#' @param .seq.col Name of the column with clonotypes' CDR3 sequences.
#' 
#' @return Data frame or a list with data frames with updated counts and proportion columns
#' and rows with unique clonotypes only.
#' 
#' @examples
#' \dontrun{
#' tmp <- data.frame(A = c('a','a','b','c', 'a')
#' B = c('V1', 'V1','V1','V2', 'V3')
#' C = c(10,20,30,40,50), stringsAsFactors = F)
#' tmp
#' #   A  B  C
#' # 1 a V1 10
#' # 2 a V1 20
#' # 3 b V1 30
#' # 4 c V2 40
#' # 5 a V3 50
#' group.clonotypes(tmp, 'B', 'C', 'A')
#' #  A  B  C
#' #  1 a V1 30
#' #  3 b V1 50
#' #  4 c V2 30
#' #  5 a V3 40
#' group.clonotypes(tmp, NA, 'C', 'A')
#' #   A  B  C
#' # 1 a V1 80
#' # 3 b V1 30
#' # 4 c V2 40
#' # For tcR data frame:
#' data(twb)
#' twb1.gr <- group.clonotypes(twb[[1]])
#' twb.gr <- group.clonotypes(twb)
#' }
group.clonotypes <- function (.data, .gene.col = 'V.gene', .count.col = 'Read.count',
                              .prop.col = 'Read.proportion', .seq.col = 'CDR3.amino.acid.sequence') {
  if (has.class(.data, 'list')) {
    return(lapply(.data, group.clonotypes, .gene.col = .gene.col, .count.col = .count.col, .seq.col = .seq.col))
  }
  
  namesvec <- c(.seq.col)
  if (!is.na(.gene.col)) {
    namesvec <- c(namesvec, .gene.col)
  }
  val <- dplyr::summarise_(dplyr::grouped_df(.data, lapply(namesvec, as.name)), value = paste0('sum(', .count.col, ')', sep = '', collapse = ''))$value
  .data <- .data[!duplicated(.data[, namesvec]), ]
  .data[, .count.col] <- val
  .data[, .prop.col] <- val / sum(val)
  .data
}


#' Shuffling data frames.
#' 
#' @aliases permutedf unpermutedf
#' 
#' @description
#' Shuffle the given data.frame and order it by the Read.count column or un-shuffle
#' a data frame and return it to the initial order.
#' 
#' @usage
#' permutedf(.data)
#' 
#' unpermutedf(.data)
#' 
#' @param .data MiTCR data.frame or list of such data frames.
#' 
#' @return Shuffled data.frame or un-shuffled data frame if \code{.data} is a data frame, else list of such data frames.
permutedf <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, permutedf))
  }
  shuffle<-.data[sample(nrow(.data)),]
  shuffle[order(shuffle$Read.count, decreasing=T),]
}

unpermutedf <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, unpermutedf))
  }
  .data[do.call(order, .data),]
}


#' Get a random subset from a data.frame.
#' 
#' @description
#' Sample rows of the given data frame with replacement.
#' 
#' @param .data Data.frame or a list with data.frames
#' @param .n Sample size if integer. If in bounds [0;1] than percent of rows to extract. "1" is a percent, not one row!
#' @param .replace if T then choose with replacement, else without.
#' 
#' @return Data.frame of nrow .n or a list with such data.frames.
sample.clones <- function (.data, .n, .replace = T) { 
  if (has.class(.data, 'list')) { return(lapply(.data, sample.clones, .n = .n)) }
  .data[sample(1:nrow(.data), if (.n > 1) .n else round(nrow(.data) * .n), replace = .replace), ]
}


#' Check if a given object has a given class.
#' 
#' @param .data Object.
#' @param .class String naming a class.
#' 
#' @return Logical.
has.class <- function (.data, .class) { .class %in% class(.data) }


#' Copy the up-triangle matrix values to low-triangle.
#' 
#' @param mat Given up-triangle matrix.
#' 
#' @return Full matrix.
matrixdiagcopy<-function(mat){
  for (i in 1:ncol(mat))
    for (j in 1:nrow(mat))
      if (i<j)
      {
        mat[j,i]<-mat[i,j]
      }
  mat
}


#' Simple functions for manipulating progress bars.
#' 
#' @aliases add.pb set.pb
#' 
#' @description
#' Set the progress bar with the given length (from zero) or add value to the created progress bar.
#' 
#' @usage
#' set.pb(.max)
#' 
#' add.pb(.pb, .value = 1)
#' 
#' @param .max Length of the progress bar.
#' @param .pb Progress bar object.
#' @param .value Value to add to the progress bar.
#' 
#' @return Progress bar (for set.pb) or length-one numeric vector giving the previous value (for add.pb).
set.pb <- function (.max) {
  txtProgressBar(min=0, max=.max, style=3)
}

add.pb <- function (.pb, .value = 1) {
  setTxtProgressBar(.pb, .pb$getVal() + .value)
}


#' Get a sample from matrix with probabilities.
#' 
#' @description
#' Get a sample from matrix or data frame with pair-wise probabilities.
#' 
#' @param .table Numeric matrix or data frame with probabilities and columns and rows names.
#' @param .count Number of sample to fetch.
#' 
#' @return Character matrix with nrow == .count and 2 columns. row[1] in row.names(.table), row[2] in colnames(.table).
sample2D <- function (.table, .count = 1) {
  if (has.class(.table, 'data.frame')) {
    .table <- as.matrix(.table)
  }
  melted.table <- melt(.table)
  as.matrix(melted.table[sample(1:nrow(melted.table), .count, T, melted.table[,3]),c(1,2)])
}


#' Apply function to every pair of data frames from a list.
#' 
#' @aliases apply.symm apply.asymm
#' 
#' @description
#' Apply the given function to every pair in the given datalist. Function either
#' symmetrical (i.e. fun(x,y) == fun(y,x)) or assymmetrical (i.e. fun(x,y) != fun(y,x)).
#' 
#' @usage 
#' apply.symm(.datalist, .fun, ..., .diag = NA, .verbose = T)
#' 
#' apply.asymm(.datalist, .fun, ..., .diag = NA, .verbose = T)
#' 
#' @param .datalist List with some data.frames.
#' @param .fun Function to apply, which return basic class value.
#' @param ... Arguments passsed to .fun.
#' @param .diag Either NA for NA or something else != NULL for .fun(x,x).
#' @param .verbose if T then output a progress bar.
#' 
#' @return Matrix with values M[i,j] = fun(datalist[i], datalist[j])
#' 
#' @examples
#' \dontrun{
#' # equivalent to intersectClonesets(immdata, 'a0e')
#' apply.symm(immdata, intersectClonesets, .type = 'a0e')
#' }
apply.symm <- function (.datalist, .fun, ..., .diag = NA, .verbose = T) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set.pb(length(.datalist)^2 / 2 + length(.datalist)/2)
  for (i in 1:length(.datalist)) 
    for (j in i:length(.datalist)) {
      if (i == j && is.na(.diag)) { res[i,j] <- NA }
      else { res[i,j] <- .fun(.datalist[[i]], .datalist[[j]], ...) }
      if (.verbose) add.pb(pb)
    }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  matrixdiagcopy(res)
}

apply.asymm <- function (.datalist, .fun, ..., .diag = NA, .verbose = T) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set.pb(length(.datalist)^2)
  for (i in 1:length(.datalist)) 
    for (j in 1:length(.datalist)) {
      if (i == j && is.na(.diag)) { res[i,j] <- NA }
      else { res[i,j] <- .fun(.datalist[[i]], .datalist[[j]], ...) }
      if (.verbose) add.pb(pb)
    }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  res
}


#' Check for adequaty of distrubution.
#' 
#' @description
#' Check if the given .data is a distribution and normalise it if necessary with optional laplace correction.
#' 
#' @param .data Numeric vector of values.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. if T then do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for laplace correction.
#' @param .na.val Replace all NAs with this value.
#' @param .warn.zero if T then the function checks if in the resulted vector (after normalisation)
#' are any zeros, and print a warning message if there are some.
#' @param .warn.sum if T then the function checks if the sum of resulted vector (after normalisation)
#' is equal to one, and print a warning message if not.
#' 
#' @return Numeric vector.
check.distribution <- function (.data, .do.norm = NA, .laplace = 1, .na.val = 0, .warn.zero = F, .warn.sum = T) {
  if (sum(is.na(.data)) == length(.data)) {
    warning("Error! Input vector is completely filled with NAs. Check your input data to avoid this. Returning vectors with zeros.\n")
    return(rep.int(0, length(.data)))
  }
  
  if (is.na(.do.norm)) {
    .data[is.na(.data)] <- .na.val
    if (sum(.data) != 1) {
      .data <- .data + .laplace
      .data <- prop.table(.data + .laplace)
    }
  } else if (.do.norm) {
    .data[is.na(.data)] <- .na.val
    .data <- prop.table(.data + .laplace)
  }
  
  if (.warn.zero && (0 %in% .data)) {
	  warningText <- paste("Warning! There are", sum(which(.data == 0)), "zeros in the input vector. Function may produce incorrect results.\n")
	  if(.laplace != 1){
		  warningText <- paste(warningText, "To fix this try to set .laplace = 1 or any other small number in the function's parameters\n")
	  }else{
		  warningText <- paste(warningText, "To fix this try to set .laplace to any other small number in the function's parameters\n")
	  } 
	  warning(warningText)
  }
  
  if (.warn.sum && sum(.data) != 1) {
	  warningText <- "Warning! Sum of the input vector is NOT equal to 1. Function may produce incorrect results.\n"
	  if(!isTRUE(.do.norm)) warningText <- paste(warningText, "To fix this try to set .do.norm = TRUE in the function's parameters.\n")
	  warning(warningText)
    if (abs(sum(.data) - 1) < 1e-14) {
      message("Note: difference between the sum of the input vector and 1 is ", (sum(.data) - 1), ", which may be caused by internal R subroutines and may not affect the result at all.\n")
    }
  }
  
  .data
}


#' Get samples from a repertoire slice-by-slice or top-by-top and apply function to them.
#' 
#' @aliases top.fun slice.fun
#' 
#' @description
#' Functions for getting samples from data frames either by consequently applying
#' head functions (\code{top.fun}) or by getting equal number of rows in the moving window (\code{slice.fun})
#' and applying specified function to this samples.
#' 
#' @usage
#' top.fun(.data, .n, .fun, ..., .simplify = T)
#' 
#' slice.fun(.data, .size, .n, .fun, ..., .simplify = T)
#' 
#' @param .data Data.frame, matrix, vector or any enumerated type or a list of this types.
#' @param .n Vector of values passed to head function for top.fun or the number of slices for slice.fun.
#' @param .size Size of the slice for sampling for slice.fun.
#' @param .fun Funtions to apply to every sample subset. First input argument is a data.frame, others are passed as \code{...}.
#' @param ... Additional parameters passed to the .fun.
#' @param .simplify if T then try to simplify result to a vector or to a matrix if .data is a list.
#' 
#' @return List of length length(.n) for top.fun or .n for slice.fun.
#' 
#' @examples
#' \dontrun{
#' # Get entropy of V-usage for the first 1000, 2000, 3000, ... clones.
#' res <- top.fun(immdata[[1]], 1000, entropy.seg)
#' # Get entropy of V-usage for the interval of clones with indices [1,1000], [1001,2000], ...
#' res <- top.fun(immdata[[1]], 1000, entropy.seg)
#' }
top.fun <- function (.data, .n, .fun, ..., .simplify = T) {
  if (has.class(.data, 'list')) { 
    res <- lapply(.data, function (d) top.fun(d, .n, .fun, ..., .simplify = .simplify))
    if (.simplify) {
      res <- as.matrix(data.frame(res))
    }
    return(res)
  }
  
  res <- lapply(.n, function (nval) .fun(head(.data, nval), ...))
  names(res) <- .n
  if (.simplify) { res <- do.call(c, res) }
  res
}

slice.fun <- function(.data, .size, .n, .fun, ..., .simplify = T) {
  if (has.class(.data, 'list')) {
    res <- lapply(.data, function(d) { slice.fun(d, .size, .n, .fun, ..., .simplify = .simplify) })
    if (.simplify) {
      res <- as.matrix(data.frame(res))
    }
    return(res)
  }
  
  res <- lapply(1:.n, function (i) {
    i <- i - 1
    down <- i * .size + 1
    up <- (i + 1) * .size
    .fun(.data[down:up, ], ...)
  })
  names(res) <- sapply(1:.n, function (i) {
    paste0((i - 1) * .size + 1, ':', i * .size)
  })
  if (.simplify) { res <- do.call(c, res) }
  res
}


#' Internal function. Add legend to a grid of plots and remove legend from all plots of a grid.
#' 
#' @description
#' Given a list of ggplot2 plots, remove legend from each of them and return 
#' grid of such plots plus legend from the first vis. Suitable for plots
#' with similar legends.
.add.legend <- function (.vis.list, .vis.nrow = 2, .legend.ncol = 1) {
  leg <- gtable_filter(ggplot_gtable(ggplot_build(.vis.list[[1]] + guides(fill=guide_legend(ncol=.legend.ncol)))), "guide-box")
  grid.arrange(do.call(arrangeGrob, c(.vis.list, nrow = .vis.nrow)), leg, widths=unit.c(unit(1, "npc") - leg$width, leg$width), nrow = 1, top ='Top crosses')
}


#' Get all values from the matrix corresponding to specific groups.
#' 
#' @description 
#' Split all matrix values to groups and return them as a data frame with two columns: for values and for group names.
#' 
#' @param .mat Input matrix with row and columns names.
#' @param .groups Named list with character vectors for names of elements for each group.
#' @param .symm If T than remove symmetrical values from the input matrix.
#' @param .diag If .symm if T and .diag is F than remove diagonal values.
#' 
#' @seealso \link{repOverlap}, \link{vis.group.boxplot}
#' 
#' @examples 
#' \dontrun{
#' data(twb)
#' ov <- repOverlap(twb)
#' sb <- matrixSubgroups(ov, list(tw1 = c('Subj.A', 'Subj.B'), tw2 = c('Subj.C', 'Subj.D')));
#' vis.group.boxplot(sb)
#' }
matrixSubgroups <- function (.mat, .groups = NA, .symm = T, .diag = F) {
  
  .intergroup.name <- function (.gr1, .gr2) {
    tmp <- sort(c(.gr1, .gr2))
    paste0(tmp[1], ':', tmp[2])
  }
  
  if (.symm) {
    .mat[lower.tri(.mat, !.diag)] <- NA
  }
  
  data <- melt(.mat, na.rm = T)
  names(data) <- c("Rep1", "Rep2", "Value")
  data$Group <- 'no-group'
  if (!is.na(.groups[1])) {
    for (i in 1:length(.groups)) {
      for (j in 1:length(.groups)) {
        gr1 <- .groups[[i]]
        gr2 <- .groups[[j]]
        gr1.name <- names(.groups)[i]
        gr2.name <- names(.groups)[j]
        if (gr1.name == gr2.name) {
          data$Group[data$Rep1 %in% gr1 & data$Rep2 %in% gr2] <-
            gr1.name
        } else {
          if (!(.intergroup.name(gr2.name, gr1.name) %in% data$Group)) {
            data$Group[data$Rep1 %in% gr1 & data$Rep2 %in% gr2] <-
              .intergroup.name(gr1.name, gr2.name)
            data$Group[data$Rep2 %in% gr1 & data$Rep1 %in% gr2] <-
              .intergroup.name(gr1.name, gr2.name)
          }
        }
      }
    }
  }
  
  data[, c('Group', 'Value')]
}


#' Compute the Euclidean distance among principal components.
#' 
#' @description 
#' Compute the Euclidean distance among principal components.
#' 
#' @param .pcaobj An object returned by \code{prcomp}.
#' @param .num.comps On how many principal components compute the distance.
#' 
#' @return Matrix of distances.
#' 
#' @seealso \link{prcomp}, \link{pca.segments}, \link{repOverlap}, \link{permutDistTest}
#' 
#' @examples
#' \dontrun{
#' mat.ov <- repOverlap(AS_DATA, .norm = T)
#' mat.gen.pca <- pca.segments(AS_DATA, T, .genes = HUMAN_TRBV)
#' mat.ov.pca <- prcomp(mat.ov, scale. = T)
#' mat.gen.pca.dist <- pca2euclid(mat.gen.pca)
#' mat.ov.pca.dist <- pca2euclid(mat.ov.pca)
#' permutDistTest(mat.gen.pca.dist, list(<list of groups here>))
#' permutDistTest(mat.ov.pca.dist, list(<list of groups here>))
#' }
pca2euclid <- function (.pcaobj, .num.comps = 2) {
  mat <- .pcaobj$x
  if (.num.comps > 0) { mat <- mat[,1:.num.comps] }
  as.matrix(dist(mat))
}