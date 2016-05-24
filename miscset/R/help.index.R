#' @name help.index
#' @keywords help index browser
#' @author Sven E. Templer
#' @title Open The Package Help Index Page in a Browser
#' @description 
#' Given a package name or string, start the package help index page
#' in a browser.
#' @param pkg A character string or expression with the name of a package.
#' @param browser The browser to display. \code{text} and \code{pdf} don't
#' use a browser, but builtin text/pdf (help_type). Otherwise a character
#' string for the browser program binary to call or function.

#' @export help.index
help.index <- function (pkg, browser = NULL) {
  
  pkg <- as.character(substitute(pkg))
  ht <- getOption("help_type")
  
  if (!is.null(browser)) {
    if (is.function(browser)) {
      ht <- "html"
      options(browser=browser)
    } else if (is.character(browser)) {
      ht <- switch(browser, text=, pdf=browser, "html")
      browser <- switch(
        browser,
        html = getOption("browser"), 
        rstudio = function (x) .Call("rs_browseURL", url), 
        browser)
      options(browser=browser)
    }
  }
  
  help(package = (pkg), help_type = ht)
  
}
