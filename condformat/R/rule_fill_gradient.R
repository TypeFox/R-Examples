#' Fill column with sequential colour gradient
#'
#' Fills the background color of a column using a gradient based on
#' the values given by an expression
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
#' @inheritParams scales::seq_gradient_pal
#' @param limits range of limits that the gradient should cover
#' @param na.value fill color for missing values
#' @param lockcells logical value determining if no further rules should be applied to the affected cells.
#'
#' @return The condformat_tbl object, with the added formatting information
#' @examples
#' data(iris)
#' condformat(iris[c(1:5, 70:75, 120:125), ]) + rule_fill_gradient(Sepal.Length)
#' condformat(iris[c(1:5, 70:75, 120:125), ]) +
#'   rule_fill_gradient(Species, expression=Sepal.Length - Sepal.Width)
#' @export
#' @importFrom lazyeval lazy_dots lazy
rule_fill_gradient <- function(...,
                               expression,
                               low = "#132B43", high = "#56B1F7",
                               space = "Lab",
                               na.value = "#7F7F7F",
                               limits=NA,
                               lockcells=FALSE) {
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
                         low = low, high = high, space = space, na.value = na.value,
                         limits = limits, lockcells = lockcells),
                    class = c("condformat_rule", "rule_fill_gradient"))
  return(rule)
}

#' @rdname rule_fill_gradient
#'
#' @family rule
#' @param columns [SE] a character vector with the column names (Only in rule_fill_gradient_)
#' @param env [SE] the environment where `expression` is to be evaluated (Only in rule_fill_gradient_)
#' @export
#' @examples
#' data(iris)
#' condformat(iris) + rule_fill_gradient_(columns=c("Sepal.Length"))
#' condformat(iris) + rule_fill_gradient_("Species", expression="Sepal.Length-Sepal.Width")
rule_fill_gradient_ <- function(columns,
                                expression,
                                low = "#132B43", high = "#56B1F7",
                                space = "Lab",
                                na.value = "#7F7F7F",
                                limits = NA,
                                lockcells = FALSE,
                                env=parent.frame()) {
  if (missing(expression)) {
    if (length(columns) > 1) {
      warning("rule_fill_gradient_ applied to multiple variables, using the first given variable as expression")
    }
    expression <- columns[1]
  }
  rule <- structure(list(columns = columns,
                         expression = expression, env = env,
                         low = low, high = high, space = space, na.value = na.value,
                         limits = limits, lockcells = lockcells),
                    class = c("condformat_rule", "rule_fill_gradient_"))
  return(rule)
}

#' @importFrom dplyr select_vars_
#' @importFrom lazyeval lazy_eval
applyrule.rule_fill_gradient <- function(rule, finalformat, xfiltered, xview, ...) {
  columns <- dplyr::select_vars_(colnames(xview), rule$columns)
  values_determining_color <- lazyeval::lazy_eval(rule$expression, xfiltered)
  rule_fill_gradient_common(rule, finalformat, xview, columns, values_determining_color)
}

applyrule.rule_fill_gradient_ <- function(rule, finalformat, xfiltered, xview, ...) {
  columns <- rule$columns
  values_determining_color <- eval(rule$expression, envir = xfiltered, enclos = rule$env)
  rule_fill_gradient_common(rule, finalformat, xview, columns, values_determining_color)
}

#' @importFrom scales seq_gradient_pal rescale
#' @importFrom assertthat are_equal
rule_fill_gradient_common <- function(rule, finalformat, xview,
                                      columns, values_determining_color) {
  if (identical(rule$limits, NA)) {
    limits <- range(values_determining_color, na.rm = TRUE)
  } else {
    limits <- rule$limits
  }

  col_scale <- scales::seq_gradient_pal(low = rule$low, high = rule$high, space = rule$space)

  values_rescaled <- scales::rescale(x = values_determining_color, from = limits)
  colours_for_values <- col_scale(values_rescaled)
  assertthat::are_equal(length(colours_for_values), nrow(xview))
  colours_for_values <- matrix(colours_for_values,
                               nrow = nrow(xview), ncol = ncol(xview), byrow = FALSE)

  finalformat <- fill_css_field_by_cols(finalformat, "background-color",
                                        colours_for_values, columns,
                                        xview, rule$lockcells)
  return(finalformat)
}
