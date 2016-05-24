scanRepeated <-  function (file, what, ...) {
    lines <- readLines(file)
    scanText <- function(text, what, quiet = TRUE, ...) {
        con <- textConnection(text, "r")
        value <- scan(con, what, quiet = quiet, ...)
        close(con)
        value
    }
    mapply(scanText, lines, what, MoreArgs = list(...),
           SIMPLIFY = FALSE)
}
