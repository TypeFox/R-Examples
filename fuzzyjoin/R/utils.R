"%||%" <- function(x, y) if (is.null(x)) y else x


common_by <- function(by = NULL, x, y) {
  if (is.list(by)) return(by)

  if (!is.null(by)) {
    by <- by[!duplicated(by)]
    x <- names(by) %||% by
    y <- unname(by)

    # If x partially named, assume unnamed are the same in both tables
    x[x == ""] <- y[x == ""]

    return(list(x = x, y = y))
  }

  by <- intersect(tbl_vars(x), tbl_vars(y))
  if (length(by) == 0) {
    stop("No common variables. Please specify `by` param.", call. = FALSE)
  }
  message("Joining by: ", utils::capture.output(dput(by)))

  list(
    x = by,
    y = by
  )
}
