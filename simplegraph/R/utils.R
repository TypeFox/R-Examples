
merge_named_lists <- function(list1, list2, default = character()) {

  names <- unique(sort(c(names(list1), names(list2))))

  res <- lapply(names, function(n) { c(list1[[n]], list2[[n]]) %||% default })

  names(res) <- names
  res
}

`%||%` <- function(l, r) if (is.null(l)) r else l

#' Create a data frame, more robust than \code{data.frame}
#'
#' It does not create factor columns.
#' It recycles columns to match the longest column.
#'
#' @param ... Data frame columns.
#' @return The constructed data frame.
#'
#' @keywords internal

data_frame <- function(...) {

  args <- list(...)

  ## Replicate arguments if needed
  len <- vapply(args, length, numeric(1))
  stopifnot(length(setdiff(len, 1)) <= 1)
  len <- max(0, max(len))
  args <- lapply(args, function(x) rep(x, length.out = len))

  ## Names
  names <- as.character(names(args))
  length(names) <- length(args)
  names <- ifelse(
    is.na(names) | names == "",
    paste0("V", seq_along(args)),
    names)

  structure(args,
            class = "data.frame",
            names = names,
            row.names = seq_along(args[[1]]))
}

drop_rownames <- function(df) {
  row.names(df) <- NULL
  df
}
