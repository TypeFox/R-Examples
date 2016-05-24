
.onAttach <- function(...) {
    msg <- "From version 0.3 functions tanorm and glyco.outliers expect data frames in long format."
    packageStartupMessage(paste(strwrap(msg), collapse = "\n"))
    return(TRUE)
}
