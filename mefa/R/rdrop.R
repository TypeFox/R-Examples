rdrop <-
function(x, cutoff=0, attrib=FALSE)
{
    if (inherits(x, "mefa")) {
        rs <- rowSums(xtab(x))
    } else if (inherits(x, "Mefa")) {
        rs <- rowSums(Matrix::as.matrix(mefa4::xtab(x)))
    } else rs <- rowSums(x)
    if (any(rs <= cutoff)) {
        exclude <- which(rs <= cutoff)
        rval <- x[-exclude,]
        if (attrib) {
            attr(rval, "exclude") <- exclude
            attr(attr(rval, "exclude"), "cutoff") <- cutoff
            attr(attr(rval, "exclude"), "margin") <- 1
        }
        rval
    } else x
}
