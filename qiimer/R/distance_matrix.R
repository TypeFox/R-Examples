#' Read a QIIME distance matrix file.
#'
#' @param filepath Path to QIIME distance matrix file.
#' @return A distance matrix.
#' @export
read_qiime_distmat <- function(filepath) {
  distmat_file <- file(filepath, 'rt')
  header_line <- readLines(distmat_file, n=1)
  column_names <- unlist(strsplit(header_line, "\t"))

  # Provide column classes to speed up conversion.
  column_classes <- rep("numeric", times=length(column_names))
  # From the read.table docs: Note that colClasses is specified per
  # column (not per variable) and so includes the column of row names.
  column_classes[1] <- "character"
  
  distmat <- as.matrix(read.table(
    distmat_file, col.names=column_names, row.names=1,
    colClasses=column_classes, sep="\t", quote=""))

  close(distmat_file)
  as.dist(distmat)
}

#' Retrieve distances from a `"dist"` object.
#' 
#' @param d A distance matrix object of class `"dist"`.
#' @param idx1,idx2 Indicies specifying the distances to extract.
#' @return A vector of distances.
#' @export
#' @examples
#' data(relmbeta_dist)
#' dist_get(relmbeta_dist, "A1", "A2")
#' dist_get(relmbeta_dist, "A1", c("A2", "A3", "A4", "A5"))
#' dist_get(relmbeta_dist, c("A1", "A2", "A3"), c("B1", "B2", "B3"))
dist_get <- function (d, idx1, idx2) {
  d <- as.dist(d)
  if (is.character(idx1)) {
    idx1 <- match(idx1, attr(d, "Labels"))
  }
  if (is.character(idx2)) {
    idx2 <- match(idx2, attr(d, "Labels"))
  }
  n <- attr(d, "Size")
  if (any(is.na(idx1) | (idx1 < 1) | (idx1 > n))) {
    stop("idx1 out of range")
  }
  if (any(is.na(idx2) | (idx2 < 1) | (idx2 > n))) {
    stop("idx2 out of range")
  }
  i <- pmin(idx1, idx2)
  j <- pmax(idx1, idx2)
  # Zeros are eliminated from index vectors
  # Need to fill with NA if i and j are equal
  idx <- ifelse(i == j, NA, n*(i-1) - i*(i-1)/2 + j-i)
  ifelse(i == j, 0, d[idx])
}

#' Extract parts of a `"dist"` object.
#' 
#' @param d A distance matrix object of class `"dist"`.
#' @param idx Indices specifying the subset of distances to extract.
#' @return A distance matrix.
#' @export
#' @examples
#' data(relmbeta_dist)
#' dist_subset(relmbeta_dist, c("A1", "A2", "A3", "A4", "A5"))
dist_subset <- function (d, idx) {
  as.dist(as.matrix(d)[idx, idx])
}

#' Create a data frame of distances between groups of items.
#'
#' @param d A distance matrix object of class `"dist"`.
#' @param g A factor representing the groups of objects in `d`.
#' @return A data frame with 6 columns. "Item1" and "Item2" identify the
#'   items compared, using the label if available. Likewise, "Group1" and 
#'   "Group2" identify the groups of the items. "Label" is a factor giving a
#'   convenient label for the type of comparison. Finally, "Distance" contains
#'   the distance of interest.
#' @export
#' @examples
#' data(relmbeta_dist)
#' data(relmbeta)
#' head(dist_groups(relmbeta_dist, relmbeta$Diet))
dist_groups <- function(d, g) {
  d <- as.dist(d)
  g <- as.factor(g)
  dsize <- attr(d, "Size")
  if (length(g) != dsize) {
    stop(
      "Length of grouping vector (g) must equal number of observations in ",
      "dist object (d)")
  }
  dlabels <- attr(d, "Labels")
  idxs <- combn(dsize, 2)
  idx1 <- idxs[1,]
  idx2 <- idxs[2,]
  
  # For the between group labels, we need to keep the groups in factor order.
  # Here, we record the level of the group to use for the first and second 
  # parts of the label.
  level1 <- levels(g)[pmin(as.numeric(g[idx1]), as.numeric(g[idx2]))]
  level2 <- levels(g)[pmax(as.numeric(g[idx1]), as.numeric(g[idx2]))]

  data.frame(
    Item1 = if (is.null(dlabels)) idx1 else dlabels[idx1],
    Item2 = if (is.null(dlabels)) idx2 else dlabels[idx2],
    Group1 = g[idx1],
    Group2 = g[idx2],
    Label = ifelse(
      level1 == level2, 
      paste("Within", level1), 
      paste("Between", level1, "and", level2)),
    Distance = dist_get(d, idx1, idx2))
}
