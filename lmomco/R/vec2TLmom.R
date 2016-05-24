"vec2TLmom" <-
function(vec, ...) {
    z <- vec2lmom(vec, ..., checklmom=FALSE)
    z$source <- "vec2TLmoms"
    return(z)
}
