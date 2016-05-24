`mefaCompare` <-
function(x1, x2, strict = FALSE)
{
    if (!inherits(x1, "mefa") || !inherits(x2, "mefa"))
        stop("both compared objects must be of class 'mefa'")
    if (!identical(dim(x1), dim(x2)))
        return(FALSE) else {

    k <- dim(x1)[3]
    xt1 <- xt2 <- list()
    if (!strict) {
        xt1[[1]] <- x1$xtab[order(rownames(x1$xtab)), order(colnames(x1$xtab))]
        xt2[[1]] <- x2$xtab[order(rownames(x2$xtab)), order(colnames(x2$xtab))]
    } else {
        xt1[[1]] <- x1$xtab
        xt2[[1]] <- x2$xtab
    }
    sego <- match(dimnames(x1)$segm, dimnames(x2)$segm)
    if (any(is.na(sego)))
        return(FALSE) else {
        if (k > 1)
            for (i in 1:k) {
                if (!strict) {
                xt1[[(i + 1)]] <- x1$segm[[i]][order(rownames(x1$segm[[i]])),
                    order(colnames(x1$segm[[i]]))]
                xt2[[(i + 1)]] <- x2$segm[[sego[i]]][order(rownames(x2$segm[[sego[i]]])),
                    order(colnames(x2$segm[[sego[i]]]))]
                } else {
                xt1[[(i + 1)]] <- x1$segm[[i]]
                xt2[[(i + 1)]] <- x2$segm[[sego[i]]]
                }
            }
        j <- length(xt1)
        rv <- logical(3 * j + 1)
        for (i in 1:j) {
            rv[i] <- all(xt1[[i]] == xt2[[i]])
            rv[(i + j)] <- all(rownames(xt1[[i]]) == rownames(xt2[[i]]))
            rv[(i + 2*j)] <- all(colnames(xt1[[i]]) == colnames(xt2[[i]]))
        }
#this lina can be left out due to match/sego line, but double chank is the sure
        rv[(3 * j + 1)] <- all(dimnames(x1)[[3]] %in% dimnames(x2)[[3]])
        return(all(rv))
    }}
}

