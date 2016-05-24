#' Simple Formula Interface for Grouped Nonlinear Functions
#'
#' This simple formula interface handles formulae of the form
#'   \code{dependent ~ FUN(independent, parms) | group1 + group2 + ...}.
#'
#' This function is used by \code{\link{all_growthmodels}} and normally not
#'   called for the user.
#'
#' @param formula a model formula specifying dependent and
#'   independent variables, nonlinear model and grouping variables in the form:
#'   \code{dependent ~ FUN(independent, parms) | group1 + group2 + ...}.
#'   FUN can be a name of an existing growth model (e.g. \code{grow_logistic})
#'   or a valid user-defined function (see \code{\link{growthmodel}}).
#'
#' @return a list with the elements \code{FUN}, \code{valuevar},  \code{timevar},
#' and \code{groups}
#'
#' @seealso \code{\link{multisplit}},  \code{\link{split}},  \code{\link{parse_formula}}
#'
#' @examples
#'
#' ret <- parse_formula_nonlin(y ~ f(x, parms) | a + b + c)
#'
#' @keywords internal
#'
#' @export
#'
parse_formula_nonlin <- function(formula) {

  form <- as.formula(formula)

  valuevar <- as.character(form[[2]])  # dependent variable
  RHS      <- form[[3]]                # f(time, parms) | group1 + group2 + ...

  rhs <- as.character(RHS)

  if (rhs[1] == "|") {    # with grouping

    FUN <- gsub("^\\s+|\\s+$", "", rhs[[2]]) # trim

    ## example: length(grep("^.*[(].*[)]$", "test(x, y)"))

    if (length(grep("^.*[(].*[)]$", FUN)) < 1) {     # no nonlinear part
      timevar <- rhs[2]
      FUN1 <- FUN2 <- NULL
    } else {                                         # with nonlinear part

      FUN1 <- parse(text = FUN)                      # full expression
      #FUN2 <- parse(text = gsub("[(].*", "", FUN))  # function name only
      FUN2 <- gsub("[(].*", "", FUN)                 # as character

      vars <- all.vars(FUN1)
      if (length(vars) != 2)
        stop ("Nonlinear part of the formula should be in the form FUN(time, parms)")
      timevar <- vars[1]
      # vars[2] is just a dummy placeholder
    }

    groups   <- gsub("[*:]", "+", rhs[3])       # convert "*" or ":" to "+"
    groups   <- unlist(strsplit(groups, "[+]")) # split right hand side
    groups   <- gsub("^\\s+|\\s+$", "", groups) # trim

  } else {    # no grouping

    if (length(rhs) == 1) {    # rhs has only independend variable
      timevar <- rhs[1]
      FUN1 <- FUN2 <- NULL
      groups <- NULL
    } else {                   # rhs contains nonlinear function
      timevar <- rhs[2]
      FUN <- rhs[1]
      FUN1 <- parse(text = FUN)
      FUN2 <- gsub("[(].*", "", FUN)
      groups <- NULL
    }
  }

  list(FUN1 = FUN1,
       FUN2 = FUN2,
       valuevar = valuevar,
       timevar = timevar,
       groups = groups)
}
