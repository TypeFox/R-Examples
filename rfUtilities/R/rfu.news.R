#' @title rfUtilities news
#' @description Displays release notes
#' @param ...  not used
#' @export
rfu.news <- function(...) {
    newsfile <- file.path(system.file(package="rfUtilities"), "NEWS")
    file.show(newsfile)
}
