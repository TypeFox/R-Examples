#' Customizes appearance of condformat object
#'
#' @param ... Arguments to be passed to htmlTable
#' @seealso \code{\link[htmlTable]{htmlTable}}
#' @examples
#' data(iris)
#' condformat(head(iris)) + theme_htmlTable(caption="Table 1: My iris table", rnames=FALSE)
#' @export
theme_htmlTable <- function(...) {
  htmlargs <- list(...)

  theme <- structure(list(htmlargs = htmlargs),
                     class = c("theme_htmlTable", "condformat_theme"))
  return(theme)
}

render_theme.theme_htmlTable <- function(themeobj, finaltheme, xview, ...) {
  for (paramname in names(themeobj$htmlargs)) {
    finaltheme[[paramname]] <- themeobj$htmlargs[[paramname]]
  }
  finaltheme
}
