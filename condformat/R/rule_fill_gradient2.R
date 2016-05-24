#' Fill column with divergent colour gradient
#'
#' Fills the background color of a column using a three colors gradient based on
#' the values by an expression
#'
#' @family rule
#'
#' @param ... Comma separated list of unquoted column names.
#'            If expression is also given, then this list can use any of the
#'            \code{\link[dplyr]{select}} syntax possibilities.
#' @param expression An expression to be evaluated with the data frame this rule applies to.
#'                   It should evaluate to a numeric vector, that will be used
#'                   to compute the color gradient values.
#'                   In rule_fill_gradient_, a character string.
#' @inheritParams scales::div_gradient_pal
#' @param midpoint the value used for the middle color (the median by default)
#' @param limits range of limits that the gradient should cover
#' @param na.value fill color for missing values
#' @param lockcells logical value determining if no further rules should be applied to the affected cells.
#'
#' @return The condformat_tbl object, with the added formatting information
#' @examples
#' data(iris)
#' condformat(iris[c(1:5, 70:75, 120:125), ]) +
#'  rule_fill_gradient2(Sepal.Length)
#' condformat(iris[c(1:5, 70:75, 120:125), ]) +
#'  rule_fill_gradient2(Species, expression=Sepal.Length - Sepal.Width)
#' @export
#' @importFrom lazyeval lazy_dots lazy
#' @importFrom scales muted
rule_fill_gradient2 <- function(...,
                                expression,
                                low = scales::muted("red"), mid="white", high = scales::muted("blue"),
                                midpoint = NA,
                                space = "Lab",
                                na.value = "#7F7F7F",
                                limits = NA,
                                lockcells = FALSE) {
  columns <- lazyeval::lazy_dots(...)
  if (missing(expression)) {
    if (length(columns) > 1) {
      warning("rule_fill_gradient applied to multiple variables, using the first given variable as expression")
    }
    expression <- columns[[1]]
  } else {
    expression <- lazyeval::lazy(expression)
  }

  rule <- structure(list(columns = columns, expression = expression,
                         low = low, mid = mid, high = high, midpoint = midpoint, space = space,
                         na.value = na.value, limits = limits, lockcells = lockcells),
                    class = c("condformat_rule", "rule_fill_gradient2"))
  return(rule)
}

#' @rdname rule_fill_gradient2
#'
#' @family rule
#' @param columns [SE] a character vector with the column names (Only in rule_fill_gradient2_)
#' @param env [SE] the environment where `expression` is to be evaluated (Only in rule_fill_gradient2_)
#' @export
#' @examples
#' data(iris)
#' condformat(iris) + rule_fill_gradient2_(columns=c("Sepal.Length"))
#' condformat(iris) + rule_fill_gradient2_("Species", expression="Sepal.Length-Sepal.Width")
rule_fill_gradient2_ <- function(columns,
                                 expression,
                                 low = scales::muted("red"), mid="white", high = scales::muted("blue"),
                                 midpoint = NA,
                                 space = "Lab",
                                 na.value = "#7F7F7F",
                                 limits = NA,
                                 lockcells = FALSE,
                                 env=parent.frame()) {
  if (missing(expression)) {
    if (length(columns) > 1) {
      warning("rule_fill_gradient2_ applied to multiple variables, using the first given variable as expression")
    }
    expression <- columns[1]
  }
  rule <- structure(list(columns = columns,
                         expression = expression, env = env,
                         low = low, mid = mid, high = high,
                         space = space, na.value = na.value,
                         midpoint = midpoint,
                         limits = limits, lockcells = lockcells),
                    class = c("condformat_rule", "rule_fill_gradient2_"))
  return(rule)
}

#' @importFrom dplyr select_vars_
#' @importFrom lazyeval lazy_eval
applyrule.rule_fill_gradient2 <- function(rule, finalformat, xfiltered, xview, ...) {
  columns <- dplyr::select_vars_(colnames(xview), rule$columns)
  values_determining_color <- lazyeval::lazy_eval(rule$expression, xfiltered)
  rule_fill_gradient2_common(rule, finalformat, xview, columns, values_determining_color)
}

applyrule.rule_fill_gradient2_ <- function(rule, finalformat, xfiltered, xview, ...) {
  columns <- rule$columns
  values_determining_color <- eval(rule$expression, envir = xfiltered, enclos = rule$env)
  rule_fill_gradient2_common(rule, finalformat, xview, columns, values_determining_color)
}

#' @importFrom scales div_gradient_pal rescale
#' @importFrom assertthat are_equal
rule_fill_gradient2_common <- function(rule, finalformat, xview,
                                      columns, values_determining_color) {
  if (identical(rule$limits, NA)) {
    limits <- range(values_determining_color, na.rm = TRUE)
  } else {
    limits <- rule$limits
  }

  if (is.na(rule$midpoint)) {
    midpoint <- stats::median(values_determining_color, na.rm = TRUE)
  } else {
    midpoint <- rule$midpoint
  }

  col_scale <- scales::div_gradient_pal(low = rule$low, mid = rule$mid, high = rule$high, space = rule$space)

  values_rescaled <- scales::rescale_mid(x = values_determining_color,
                                         from = limits, mid = midpoint)

  colours_for_values <- col_scale(values_rescaled)
  assertthat::are_equal(length(colours_for_values), nrow(xview))
  colours_for_values <- matrix(colours_for_values,
                               nrow = nrow(xview), ncol = ncol(xview), byrow = FALSE)

  finalformat <- fill_css_field_by_cols(finalformat, "background-color",
                                        colours_for_values, columns, xview,
                                        rule$lockcells)
  return(finalformat)
}
