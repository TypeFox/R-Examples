#' Shiny app for EDGAR filing management.
#'
#' \code{RunEdgarShiny } is a shiny app which has functionality of downloading  
#' quarterly master index, edgar filings, daily index, and do sentiment analysis 
#' of 10-K statement.
#' 
#' RunEdgarShiny app has tabs: 'Get Master Index', 'Filings Info', 'Download Filings',
#' 'Get Daily EDGAR Info', and 'Sentiment Analysis of 10-K'.
#' See the directory structure at home page. App contains pop-up interactive messages, 
#' it can only be run on R-studio's Viewer pane.
#' Do not run this app on any browser.
#' 
#' @examples
#' \dontrun{
#' 
#' RunEdgarShiny()
#' }
#' 
RunEdgarShiny <- function() {
    appDir <- system.file("shiny-edgar", "edgarapp", package = "edgar")
    if (appDir == "") {
        stop("Could not find shiny-edgar directory. Try re-installing `edgar`.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")
} 
