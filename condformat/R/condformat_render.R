#' Prints the data frame in an html page and shows it.
#'
#' @param x A condformat_tbl object
#' @param ... optional arguments to print
#' @return the value returned by htmlTable
#' @examples
#' data(iris)
#' print(condformat(iris))
#' @export
print.condformat_tbl <- function(x, ...) {
  thetable <- condformat2html(x)
  print(thetable)
  invisible(x)
}


#' Converts the table to a htmlTable object
#'
#' @param x A condformat_tbl object
#' @return the htmlTable object
#' @importFrom htmlTable htmlTable
#' @examples
#' data(iris)
#' print(condformat(iris))
#' @export
condformat2html <- function(x) {
  finalshow <- render_show_condformat_tbl(x)
  xfiltered <- finalshow$xfiltered
  xview <- xfiltered[, finalshow$cols, drop = FALSE]
  rules <- attr(x, "condformat")$rules
  finalformat <- render_rules_condformat_tbl(rules, xfiltered, xview)
  # Rename the columns according to show options:
  colnames(xview) <- names(finalshow$cols)
  themes <- attr(x, "condformat")$themes
  finaltheme <- render_theme_condformat_tbl(themes, xview)
  thetable <- do.call(htmlTable::htmlTable, c(list(format(xview),
                                                   css.cell = finalformat$css_cell),
                                              finaltheme))
  return(thetable)
}


#' @importFrom knitr knit_print
#' @export
knit_print.condformat_tbl <- function(x, ...) {
  knitr::knit_print(condformat2html(x), ...)
}

render_theme_condformat_tbl <- function(themes, xview) {
  finaltheme <- list()
  for (themeobj in themes) {
    finaltheme <- render_theme(themeobj, finaltheme, xview)
  }
  return(finaltheme)
}


render_show_condformat_tbl <- function(x) {
  condformatopts <- attr(x, "condformat")

  finalshow <- list(xfiltered = x,
                    cols = colnames(x))
  names(finalshow$cols) <- colnames(x)

  # First we filter, then we select so we can
  # filter by variables not selected
  showobjs <- c(condformatopts$show$rows,
                condformatopts$show$cols)
  for (showobj in showobjs) {
    finalshow <- render_show(showobj, finalshow, finalshow$xfiltered)
  }

  return(finalshow)
}

merge_css_conditions <- function(initial_value, css_fields) {
  css_keys <- names(css_fields)
  output <- initial_value
  for (key in css_keys) {
    thisfield <- paste(key, css_fields[[key]], sep = ": ")
    output <- paste(output, thisfield, sep = "; ") # I don't care about a leading "; "
  }
  output <- matrix(output, nrow = nrow(initial_value), ncol = ncol(initial_value))
  return(output)
}

#' Renders the css matrix to format the xview table
#'
#' @param rules List of rules to be applied
#' @param xview Data frame with the rows and columns that will be printed
#' @param xfiltered Like xview, but with all the columns (rules
#'                  will use columns that won't be printed)
#' @return List with the CSS information
render_rules_condformat_tbl <- function(rules, xfiltered, xview) {

  finalformat <- list(css_fields = list(),
                      css_cell = matrix(data = "", nrow = nrow(xview), ncol = ncol(xview)),
                      css_cell_unlocked = matrix(data = TRUE,
                                                 nrow = nrow(xview),
                                                 ncol = ncol(xview)))

  for (rule in rules) {
    finalformat <- applyrule(rule, finalformat, xfiltered, xview)
  }
  if (length(finalformat$css_fields) > 0) {
    finalformat$css_cell <- merge_css_conditions(finalformat$css_cell, finalformat$css_fields)
  }
  return(finalformat)
}

render_show <- function(showobj, finalshow, x, ...) UseMethod("render_show")

render_show.default <- function(showobj, finalshow, x , ...) {
  finalshow
}

applyrule <- function(rule, finalformat, xfiltered, xview, ...) UseMethod("applyrule")


applyrule.default <- function(rule, finalformat, xfiltered, xview, ...) {
   finalformat
}

render_theme <- function(themeobj, finaltheme, xview, ...) UseMethod("render_theme")


render_theme.default <- function(themeobj, finaltheme, xview, ...) {
  finaltheme
}
