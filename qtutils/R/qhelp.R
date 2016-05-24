
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A Qt web browser widget
##' @param url The URL to view.
##' @param ... Ignored.
##' @return A QWebView instance
##' @author Deepayan Sarkar
qwebbrowser <- function(url, ...)
{
    w <- Qt$QWebView()
    w$load(Qt$QUrl(url))
    ## w$setWindowTitle(w$title())
    ## w$show()
    w
}



