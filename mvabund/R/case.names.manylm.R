case.names.manylm	<- case.names.manyglm <-
function (object, full = FALSE, ...) 
{
    w <- weights(object)
    dn <- rownames(residuals(object))
    if (full || is.null(w)) 
        dn
    else dn[w != 0]
}
