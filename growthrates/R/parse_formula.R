#' Simple Formula Interface
#'
#' This simple formula interface handles formulae of the form
#'   \code{dependent ~ independent | group1 + group2 + ...}.
#'
#' This function is used by \code{\link{multisplit}} and normally not called
#'   by the user.
#'
#' @param grouping a model formula specifying dependent,
#'   independent and grouping variables in the form:
#'   \code{dependent ~ independent | group1 + group2 + ...}.
#'
#' @return a list with the elements \code{valuevar},  \code{timevar}, and
#'   \code{groups}
#'
#' @seealso \code{\link{multisplit}},  \code{\link{split}}
#'
#' @examples
#'
#' parse_formula(y ~ x | a+b+c)
#'
#' @keywords internal
#'
#' @export
#'
parse_formula <- function(grouping) {
  tm <- terms(grouping)

  valuevar <- as.character(tm[[2]])
  RHS      <- as.character(tm[[3]])
  if (RHS[1] == "|") {
    timevar <- RHS[2]    # with grouping
  } else {
    timevar <- RHS[1]    # without grouping
  }
  groups   <-
    gsub("[*:]", "+", RHS[3])       # convert "*" or ":" to "+"
  groups   <-
    unlist(strsplit(groups, "[+]")) # split right hand side
  groups   <- gsub("^\\s+|\\s+$", "", groups) # trim

  ## replace NA by NULL
  if (is.na(groups[1])) groups <- NULL

  list(valuevar = valuevar,
       timevar = timevar,
       groups = groups)
}
