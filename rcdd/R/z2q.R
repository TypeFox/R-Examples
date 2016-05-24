z2q <- function(numer, denom, canonicalize = TRUE) {
    if (any(denom == 0))
        stop("zero denominator")
    if (length(numer) != length(denom))
        stop("length(numer) != length(denom)")
    ans <- paste(numer, denom, sep = "/")
    mostattributes(ans) <- attributes(numer)
    if (canonicalize)
        ans <- q2q(ans)
    return(ans)
}
