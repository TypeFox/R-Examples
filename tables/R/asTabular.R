as.tabular <- function(x, like = NULL) {
  UseMethod("as.tabular")
}

as.tabular.default <- function(x, like = NULL) {
  dim <- dim(x)
  if (length(dim) != 2) stop("Cannot convert to tabular")
  if (is.null(like)) {
    rows <- factor(levels=rownames(x, do.NULL = FALSE))
    cols <- factor(levels=colnames(x, do.NULL = FALSE))
    names <- names(dimnames(x))
    if (!length(names))
    names <- c("", "")

    # Make a fake table
    df <- data.frame(A=rows, B=cols)
    formula <- Factor(A, names[1]) ~ Factor(B, names[2])
    if (names[1] == "") formula[[2]] <- call("*", quote(Heading()), formula[[2]])
    if (names[2] == "") formula[[3]] <- call("*", quote(Heading()), formula[[3]])
    like <- tabular(formula, data=df)
  } else {
    if (!inherits(like, "tabular"))
      stop("'like' must be a 'tabular' object")
    if (any(dim != dim(like)))
      stop("Dimension of 'x' must match dimension of 'like'")
  }
  
  # Make x into a list, and add those attributes            
  mode(x) <- "list"
  attributes(x) <- attributes(like)
  x
}

as.tabular.data.frame <- function(x, like = NULL) {
  dimnames <- dimnames(x)
  res <- matrix(1, nrow(x), ncol(x))
  mode(res) <- "list"
  for (j in seq_len(ncol(x))) {
    col <- x[,j]
    mode(col) <- "list"
    res[,j] <- col
  }
  dimnames(res) <- dimnames
  as.tabular.default(res, like = like)
}
    