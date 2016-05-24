##
##  p o l y 2 s t r . R  Print Polynomial
##


poly2str <- function(p, svar = "x", smul = "*",
                     d = options("digits")$digits) {
    if (length(p) == 0) return("")
    if (!is.numeric(p))
        stop("Argument 'p' must be a numeric vector.")

    while (p[1] == 0 && length(p) > 1)
        p <- p[2:length(p)]
    if (length(p) == 1) return(as.character(p))

    s <- sign(p)
    p <- abs(p)

    p <- formatC(p, digits = d)
    p <- sub("^\\s+", "", p)

    n <- length(p) - 1
    S <- ""

    s1 <- if (s[1] == 1) "" else "-"
    S <- paste(s1, p[1], smul, svar, "^", n, sep = "")

    for (i in 2:(n+1)) {
        if (s[i] == 1) s1 <- " + "
        else if (s[i] == -1) s1 <- " - "
        else next

        if (n-i+1 > 1) {
            S <- paste(S, s1, p[i], smul, svar, "^", n-i+1, sep="")
        } else if (i == n) {
            S <- paste(S, s1, p[i], smul, svar, sep="")
        } else {
            S <- paste(S, s1, p[i], sep="")
        }
    }
    return(S)
}
