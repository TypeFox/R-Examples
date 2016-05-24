#' Parse a QIIME OTU table file in "calssic" format.
#'
#' @param filepath Path to OTU table file.
#' @param commented TRUE if the header line is preceeded by an additional
#'   comment line, otherwise FALSE.  This is usually the case for OTU
#'   tables generated with QIIME, so we default to TRUE.
#' @param metadata TRUE if the OTU table contains a metadata column, otherwise
#'   FALSE.  The metadata column usually contains taxonomic assignments, and
#'   must be located on the right-hand side of the table.
#' @return A list with four attributes: sample_ids, otu_ids, counts, and 
#'   metadata, a data structure similar to that returned by the python 
#'   function `qiime.parse.parse_otu_table`.  The sample_ids, otu_ids, and
#'   metadata attributes are character vectors.  The counts attribute is a
#'   matrix with one column per sample_id and one row per otu_id.
#' @export
read_qiime_otu_table <- function(filepath, commented=TRUE, metadata=TRUE) {
  f <- file(filepath, "rt")
  header_line <- readLines(f, n=1)
  if (commented) {
    header_line <- readLines(f, n=1)
  }
  col_names <- strsplit(header_line, "\t")[[1]]

  col_classes <- rep("numeric", times=length(col_names))
  col_classes[1] <- "character"
  if (metadata) {
    col_classes[length(col_classes)] <- "character"
  }

  full_otu_table <- read.table(
    f, col.names=col_names, colClasses=col_classes, sep="\t", 
    quote="", as.is=TRUE, header=FALSE)
  close(f)

  data_cols <- if (metadata) {
    2:(length(col_names) - 1) 
  } else {
    2:length(col_names)
  } 

  sample_ids <- col_names[data_cols]
  otu_ids <- as.character(full_otu_table[,1])

  counts <- as.matrix(full_otu_table[,data_cols])
  rownames(counts) <- otu_ids

  if (metadata) {
    metadata_vals <- as.character(full_otu_table[,length(col_names)])
    names(metadata_vals) <- otu_ids
  } else {
    metadata_vals <- NULL
  }
    
  list(
    sample_ids=sample_ids, otu_ids=otu_ids, counts=counts,
    metadata=metadata_vals)
}

#' Standard taxonomic ranks.
#'
#' @export
taxonomic_ranks <- c(
  "Kingdom", "Phylum", "Class", "Order", 
  "Family", "Genus", "Species")

#' Split taxonomic assignment strings
#'
#' @param assignments Character vector of taxonomic assignments.
#' @param ranks Character vector of taxonomic ranks, used as column names in the
#'   resultant data frame.
#' @param split Pattern on which to split taxa in assignment strings.
#' @param ... Additional parameters are passed to the \code{strsplit} function.
#' @return A data frame of taxonomic assignments.
#' @seealso \code{\link{taxonomic_ranks}}
#' @export
#' @examples
#' data(relmbeta_assignments)
#' a <- split_assignments(relmbeta_assignments)
#' head(a)
split_assignments <- function(assignments, ranks=taxonomic_ranks, 
  split="; ", ...) {
  a <- strsplit(as.character(assignments), split, ...)
  max_ranks <- max(sapply(a, length))
  a <- lapply(a, function (x) {
    fill_length <- max_ranks - length(x)
    c(x, rep(NA, fill_length))
  })
  a <- as.data.frame(do.call(rbind, a))
  colnames(a) <- ranks[1:ncol(a)]
  if (!is.null(names(assignments))) {
    rownames(a) <- names(assignments)
  }
  a
}

#' Reformat taxonomic assignments for presentation.
#'
#' @param assignments_df A data frame of taxonomic assignments.
#' @param rank1 The rank of taxonomy to use as the first word in the label.
#' @param rank2 The rank of taxonomy to use as the second word in the label.
#' @return A character vector of reformatted assignment labels.
#' @seealso \code{\link{split_assignments}}
#' @export
#' @examples
#' data(relmbeta_assignments)
#' a <- split_assignments(relmbeta_assignments)
#' head(simplify_assignments(a))
simplify_assignments <- function(assignments_df, rank1="Phylum", 
  rank2="Genus") {
  if (is.character(rank1)) {
    rank1 <- match(rank1, colnames(assignments_df))
  }
  if (is.character(rank2)) {
    rank2 <- match(rank2, colnames(assignments_df))
  }
  apply(assignments_df, 1, function (x) {
    x <- na.omit(as.character(x))
    n <- length(x)
    if (n == 1)     return(x)
    if (n < rank1)  return(paste(x, collapse=" "))
    if (n == rank1) return(x[rank1])
    if (n < rank2)  return(paste(x[rank1], x[length(x)]))
    return(paste(x[rank1], x[rank2]))
  })
}

#' Create a heatmap of OTU counts.
#'
#' @param otu_counts A matrix of OTU counts, one row per OTU and one column 
#'   per sample.
#' @param assignments A character vector of OTU assignments.  Length should
#'   match number of rows in otu_counts.
#' @param threshold Minimum number of OTU counts necessary for an assignment to
#'   be included in the heatmap.  Assignments are filtered after calculating
#'   the proportions, so the threshold setting does not affect the display of
#'   the remaining OTUs.
#' @param plot If true, display a plot.  If false, just return the computed
#'   abundances.
#' @param color Vector of colors to use in the heatmap.
#' @param breaks Vector of color breaks, one element greater in length than
#'   `colors`.
#' @param ... Additional arguments are passed to the pheatmap function.
#' @return A heatmap plot of the proportions of assignments in each sample,
#'   and invisibly returns a matrix of the proportions in the plot.
#' @seealso \code{\link{saturated_rainbow}}
#' @export
#' @examples
#' data(relmbeta_assignments)
#' data(relmbeta_counts)
#' a <- simplify_assignments(split_assignments(relmbeta_assignments))
#' 
#' \dontrun{
#' otu_heatmap(relmbeta_counts, a, threshold=10)
#' otu_heatmap(
#'   relmbeta_counts, a, threshold=10, 
#'   cluster_rows=FALSE, cluster_cols=FALSE, 
#'   cellwidth=12, cellheight=12)
#' }
#' 
#' heatmap_data <- otu_heatmap(relmbeta_counts, a, threshold=10, plot=FALSE)
#' head(heatmap_data)
otu_heatmap <- function(otu_counts, assignments, threshold=0, plot=T,
  color=saturated_rainbow(100),
  breaks=c(0, 1e-10, seq(0.001, 1, length.out = 100)), ...) {
  # rowsum() does not play well with factors
  assignments <- as.character(assignments)
  # NA values in assignments work fine, but produce an unnecessary warning
  assignment_counts <- suppressWarnings(rowsum(otu_counts, assignments))
  rows_to_keep <- apply(assignment_counts, 1, sum) >= threshold
  assignment_fracs <- apply(assignment_counts, 2, function (x) {x / sum(x)}) 
  assignment_fracs <- as.matrix(assignment_fracs[rows_to_keep,])
  if (plot) {
    pheatmap(assignment_fracs, color=color, breaks=breaks, ...)
  }
  # Return underlying data
  invisible(assignment_fracs)
}

#' Saturated rainbow palette.
#' 
#' This palette is specially designed for data consisting of counts.  It is
#' intended to show both presence/absence and relative proportion in the same
#' plot.  For data containing N counts in the largest sample, the saturated 
#' rainbow palette should be created with length N + 1.
#' 
#' The first element of the palette is white, indicating zero counts.  The 
#' second element is dark blue, indicating one or very few counts.  As the 
#' proportion increases within a sample, the palette transitions from 
#' blue to green, yellow, orange, and finally red.
#' 
#' The function defines a saturation limit, above which the color remains 
#' bright red.  The saturation limit is set to 40% by default, to highlight
#' items with the largest relative proportion in a sample.  The default value 
#' seems to work well for a wide range of circumstances -- it allows items that
#' are strongly dominant in a sample to be identified across the plot.  Ideally,
#' the total number of red squares should be kept low, never more than one per 
#' sample. 
#'
#' @param n Length of the palette
#' @param saturation_limit The fraction of the total palette length over which
#'   the rainbow extends.  Above this limit, the color will remain the same.
#' @return A vector of colors.
#' @export
#' @examples
#' saturated_rainbow(5)
saturated_rainbow <- function (n, saturation_limit=0.4) {
  saturated_len <- floor(n * (1 - saturation_limit))
  rainbow_colors <- rev(rainbow(n - saturated_len, start=0, end=0.6))
  last_color <- tail(rainbow_colors, n=1)
  saturated_colors <- rep(last_color, saturated_len)
  colors <- c(rainbow_colors, saturated_colors)
  colors[1] <- "#FFFFFFFF"
  colors
}
