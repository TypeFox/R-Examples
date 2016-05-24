dropFormula <- function(original, drop) {
    factors <- attr(terms(as.formula(original)),
                    "factors")
    row <-  match(drop, rownames(factors))
    whichTerms <- factors[row,] == 1
    labels <- colnames(factors)[whichTerms]
    text <- paste("~ . ",
                  paste("-", labels, collapse = " "))
    eval(parse(text = text))
}
