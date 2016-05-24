#' Fill column with discrete colors
#'
#' Fills a column or columns of a data frame using a discrete
#' colour palette, based on an expression.
#'
#' @family rule
#' @param ... Comma separated list of unquoted column names.
#'            If expression is also given, then this list can use any of the
#'            \code{\link[dplyr]{select}} syntax possibilities.
#' @param expression an expression to be evaluated with the data.
#'                   It should evaluate to a logical or an integer vector,
#'                   that will be used to determine which cells are to be coloured
#'                   In rule_fill_discrete_, a character string with the expression.
#' @param colours a character vector with colours as values and the expression
#'                possible results as names.
#' @inheritParams scales::hue_pal
#' @param na.value a character string with the CSS color to be used in missing values
#' @param lockcells logical value determining if no further rules should be applied to the affected cells.
#'
#' @return The condformat_tbl object, with the added formatting information
#' @examples
#' data(iris)
#' condformat(iris[c(1:5, 70:75, 120:125), ]) +
#'  rule_fill_discrete(Species, colours = c("setosa" = "red",
#'                                          "versicolor" = "blue",
#'                                          "virginica" = "green"))
#' condformat(iris[c(1:5, 70:75, 120:125), ]) +
#'  rule_fill_discrete(Species, expression=Sepal.Length > 4.6,
#'                     colours=c("TRUE"="red"))
#' @export
#' @importFrom lazyeval lazy_dots lazy
rule_fill_discrete <- function(...,
                               expression,
                               colours,
                               na.value = "blank",
                               h = c(0, 360) + 15, c = 100, l = 65,
                               h.start = 0, direction = 1,
                               lockcells=FALSE) {
  columns <- lazyeval::lazy_dots(...)
  if (missing(expression)) {
    if (length(columns) > 1) {
      warning("rule_fill_discrete applied to multiple variables, using the first given variable as expression")
    }
    expression <- columns[[1]]
  } else {
    expression <- lazyeval::lazy(expression)
  }

  if (missing(colours)) {
    colours = NA
  }
  rule <- structure(list(columns = columns, expression = expression,
                         colours = colours,
                         h = h, c = c, l = l, h.start = h.start, direction = direction,
                         na.value = na.value, lockcells = lockcells),
                    class = c("condformat_rule", "rule_fill_discrete"))
  return(rule)
}

#' @rdname rule_fill_discrete
#'
#' @family rule
#' @param columns [SE] a character vector with the column names (Only in rule_fill_discrete_)
#' @param env [SE] the environment where `expression` is to be evaluated (Only in rule_fill_discrete_)
#' @export
#' @examples
#' data(iris)
#' condformat(iris) + rule_fill_discrete_(columns=c("Species"))
#' condformat(iris) + rule_fill_discrete_("Species", expression="Sepal.Length > 4.6")
rule_fill_discrete_ <- function(columns,
                                expression,
                                colours,
                                h = c(0, 360) + 15, c = 100, l = 65,
                                h.start = 0, direction = 1, na.value = "blank",
                                lockcells=FALSE,
                                env=parent.frame()) {
  if (missing(expression)) {
    if (length(columns) > 1) {
      warning("rule_fill_discrete_ applied to multiple variables, using the first given variable as expression")
    }
    expression <- columns[1]
  }

  if (missing(colours)) {
    colours <- NA
  }
  rule <- structure(list(columns = columns,
                         expression = expression, env = env,
                         colours = colours,
                         h = h, c = c, l = l, h.start = h.start, direction = direction,
                         na.value = na.value, lockcells = lockcells),
                    class = c("condformat_rule", "rule_fill_discrete_"))
  return(rule)
}

#' @importFrom lazyeval lazy_eval
#' @importFrom dplyr select_vars_
applyrule.rule_fill_discrete <- function(rule, finalformat, xfiltered, xview, ...) {
  columns <- dplyr::select_vars_(colnames(xview), rule$columns)
  values_determining_color <- as.factor(lazyeval::lazy_eval(rule$expression, xfiltered))
  rule_fill_discrete_common(rule, finalformat, xfiltered, xview, columns,
                            values_determining_color)
}

applyrule.rule_fill_discrete_ <- function(rule, finalformat, xfiltered, xview, ...) {
  columns <- rule$columns
  values_determining_color <- as.factor(eval(rule$expression, envir = xfiltered, enclos = rule$env))
  rule_fill_discrete_common(rule, finalformat, xfiltered, xview, columns,
                            values_determining_color)
}

#' @importFrom scales hue_pal
#' @importFrom assertthat are_equal
rule_fill_discrete_common <- function(rule, finalformat, xfiltered, xview,
                                      columns, values_determining_color) {
  colours_for_values <- NA
  if (identical(rule$colours, NA)) {
    number_colours <- length(levels(values_determining_color))
    col_scale <- scales::hue_pal(h = rule$h, c = rule$c, l = rule$l,
                                 h.start = rule$h.start,
                                 direction = rule$direction)(number_colours)
    colours_for_values <- col_scale[as.integer(values_determining_color)]
    assertthat::are_equal(length(colours_for_values), nrow(xview))
  } else {
    colours_for_values <- rule$colours[match(values_determining_color, names(rule$colours))]
  }
  colours_for_values[is.na(colours_for_values)] <- rule$na.value
  colours_for_values <- matrix(colours_for_values,
                               nrow = nrow(xview), ncol = ncol(xview), byrow = FALSE)

  finalformat <- fill_css_field_by_cols(finalformat,
                                        "background-color", colours_for_values,
                                        columns, xview, rule$lockcells)
  return(finalformat)
}

