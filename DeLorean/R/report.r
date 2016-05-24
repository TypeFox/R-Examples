#' The filename of the R markdown report.
#'
#' @param report.name The report name
#'
#' @export
report.file <- function(report.name) {
    system.file("Rmd", sprintf("%s.Rmd", report.name),
                package="DeLorean")
}


#' The filename of the R markdown stylesheet
#'
de.lorean.stylesheet <- function() {
    system.file("Rmd", "foghorn.css", package="DeLorean")
}


#' Knit a report, the file inst/Rmd/<report.name>.Rmd must exist in
#' the package directory.
#'
#' @param dl de.lorean object
#' @param report.name The name of the report. Used to locate the R
#'  markdown report file in the package.
#'
knit.report <- function(dl, report.name) {
    knitr::knit2html(report.file(report.name),
              envir=environment(),
              stylesheet=de.lorean.stylesheet())
}


