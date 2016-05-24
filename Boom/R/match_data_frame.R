MatchDataFrame <- function(data.to.match, data.to.permute) {
  ## Given two data frames with the same data, but with rows and
  ## columns in potentially different orders, produce a pair of
  ## permutations such that data.to.permute[row.permutation, column.permutation]
  ## matches data1 as close as possible.
  ##
  ## Args:
  ##   data.to.match:  The data frame to be matched.
  ##   data.to.permute:  The data frame to be permuted.
  ##
  ## Returns:
  ##   A list containing two vectors
  ##   * column.permutation: a vector of indices such that the columns
  ##     of data.to.permute[, column.permutation] are in the same
  ##     order as the columns of data.to.match.
  ##   * row.permutation: a vector of indices such that the rows of
  ##     data.to.permute[row.permutation, ] match the rows of
  ##     data.to.match.

  if (nrow(data.to.match) != nrow(data.to.permute) ||
      ncol(data.to.match) != ncol(data.to.permute)) {
    stop("Data are of different sizes")
  }
  column.permutation <- match(colnames(data.to.match),
                              colnames(data.to.permute),
                              nomatch = NA)
  if (any(is.na(column.permutation))) {
    stop("Not all columns could be matched.")
  }

  data.to.permute <- data.to.permute[, column.permutation]
  data.to.permute.row.strings <- apply(
      data.to.permute, 1, paste0, collapse = "|")
  data.to.match.row.strings <- apply(
      data.to.match, 1, paste0, collapse = "|")
  row.permutation <- match(data.to.match.row.strings,
                           data.to.permute.row.strings,
                           nomatch = NA)
  if (any(is.na(row.permutation))) {
    stop("Not all rows could be matched.")
  }
  return(list(column.permutation = column.permutation,
              row.permutation = row.permutation))
}
