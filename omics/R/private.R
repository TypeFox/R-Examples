format.percentage <- function(x, digits) {
    paste0(format(x * 100, trim=TRUE, digits=digits, scientific=FALSE), "%")
}

response.name <- function(formula, ...) {
    tt <- terms(formula, ...)
    vars <- as.character(attr(tt, "variables"))[-1]
    vars[attr(tt, "response")]
}
