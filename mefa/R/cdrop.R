cdrop <-
function(x, cutoff=0, attrib=FALSE)
{
    if (inherits(x, "mefa")) {
        cs <- colSums(xtab(x))
    } else if (inherits(x, "Mefa")) {
        cs <- colSums(Matrix::as.matrix(mefa4::xtab(x)))
    } else cs <- colSums(x)
    if (any(cs <= cutoff)) {
        exclude <- which(cs <= cutoff)
        rval <- x[,-exclude]
        if (attrib) {
            attr(rval, "exclude") <- exclude
            attr(attr(rval, "exclude"), "cutoff") <- cutoff
            attr(attr(rval, "exclude"), "margin") <- 2
        }
        rval
    } else x
}
