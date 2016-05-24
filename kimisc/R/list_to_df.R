#' Converts a list to a name-value data frame
#'
#' This function coerces its input to a list and returns a data frame with as many
#' rows as there are list items in the input, and two columns
#' (one for the names, one for the values). If the list is not named, the
#' natural sequence will be used as item names.
#'
#' @param list_for_df The object to be converted to a data frame
#'
#' @export
list_to_df <- function(list_for_df) {
  list_for_df <- as.list(list_for_df)

  nm <- names(list_for_df)
  if (is.null(nm))
    nm <- seq_along(list_for_df)

  df <- data.frame(name = nm, stringsAsFactors = FALSE)
  df$value <- unname(list_for_df)
  df
}

#' Converts a name-value data frame to a named list
#'
#' This function converts a data frame back to a list. It is the reverse
#' operation to \link{list_to_df}.
#'
#' In a data frame with more than two columns, heuristics are applied to detect
#' the name and value column.
#'
#' @param df_for_list The data frame to be converted to a list
#'
#' @export
#' @importFrom stats setNames
df_to_list <- function(df_for_list) {
  value_cols <- which(vapply(df_for_list, is.list, logical(1L)))
  value_col <- coalesce.na(value_cols["value"], value_cols[1L])
  if (is.na(value_col)) {
    stop("No column of type list found.")
  }
  value <- df_for_list[[value_col]]

  name_cols <- setNames(nm = names(df_for_list[-value_cols]))
  name_col <- coalesce.na(name_cols["name"], name_cols[1L])

  if (!is.na(name_col)) {
    nm <- df_for_list[[name_col]]
    if (any(nm != seq_along(value)))
      names(value) <- nm
  }

  value
}
