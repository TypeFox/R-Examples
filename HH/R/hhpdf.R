## for use when running as an ordinary user
hhpdf <- function(file, ...) {invisible(NULL)}

hhdev.off <- function(...) {invisible(NULL)}

hhcapture <- function(file, text, echo=TRUE, print.eval=TRUE) {
  source(textConnection(text),
         echo=TRUE, print.eval=TRUE, keep.source=TRUE,
         max.deparse.length=500)
}

hhcode <- function(file, text) {
  cat(text)
}

hhpng <- function(file, ...) {invisible(NULL)}

hhlatex <- function(file="", ...) {
  file.tex <- Hmisc::latex(file="", ...)
  invisible(NULL)
}
