######################################################################################################################

# Function: mergeOutcomeParameter .
# Argument: two lists.
# Description: This function is used to merge two lists

mergeOutcomeParameter <- function (first, second)
{
  stopifnot(is.list(first), is.list(second))
  firstnames <- names(first)
  for (v in names(second)) {
    first[[v]] <- if (v %in% firstnames && is.list(first[[v]]) && is.list(second[[v]]))
      appendList(first[[v]], second[[v]])
    else paste0(first[[v]],' = ', lapply(second[[v]], function(x) round(x,3)))
  }
  paste0(first, collapse = ", ")
}