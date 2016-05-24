messageNAMED <- function (possiblynamed, ...) {
  if (is.null(possiblynamed)) {
    base::message(possiblynamed)
  } else {
    if (is.matrix(possiblynamed) || is.data.frame(possiblynamed)) {
      noms <- colnames(possiblynamed)
    } else noms <- names(possiblynamed)
    fieldwidth <- (max(nchar(noms), nchar(possiblynamed))) + 1
    if (!is.null(noms)) {
      fillstrings <- sapply(nchar(noms), function(n) {
        paste(rep(" ", fieldwidth - n), sep = "", collapse = "")
      })
      base::message(paste(fillstrings, noms))
    }
    if (is.data.frame(possiblynamed) || is.matrix(possiblynamed)) {
      apply(possiblynamed, 1, function(l) {
        fillstrings <- sapply(nchar(l), function(n) {
          paste(rep(" ", fieldwidth - n), sep = "", collapse = "")
        })
        base::message(paste(fillstrings, l))
      })
    } else if (is.numeric(possiblynamed)) {
      fillstrings <- sapply(nchar(possiblynamed), function(n) {
        paste(rep(" ", fieldwidth - n), sep = "", collapse = "")
      })
      base::message(paste(fillstrings, possiblynamed), ...)
    } else base::message(possiblynamed, ...)
  }
}
