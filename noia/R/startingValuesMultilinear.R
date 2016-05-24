startingValuesMultilinear <-
function (noia.multilinear, max.level = 2, max.dom = 2, e.unique = FALSE) 
{
    if (class(noia.multilinear) != "noia.multilinear") {
        stop("Object of class \"noia.multilinear\" expected\n")
    }
    a <- noia::effectsNames[2]
    d <- noia::effectsNames[3]
    e <- noia::effectsNames[4]
    ans <- list()
    nloc <- noia.multilinear$nloc
    effects <- noia.multilinear$E
    ans[["R"]] <- effects[effNames(nloc = nloc)]
    for (l1 in 1:nloc) {
        ans[[paste("a", l1, sep = "")]] <- effects[effNames(a, 
            l1, nloc)]
        if (max.dom > 0) {
            ans[[paste("d", l1, sep = "")]] <- effects[effNames(d, 
                l1, nloc)]
            if (is.na(ans[[paste("d", l1, sep = "")]])) {
                ans[[paste("d", l1, sep = "")]] <- 0
            }
        }
    }
    if ((max.level > 1) && (nloc > 1)) {
        if (e.unique) {
            ans[["ee"]] <- 0
        }
        for (l1 in 1:(nloc - 1)) {
            for (l2 in (l1 + 1):nloc) {
                if (e.unique) {
                  ans[["ee"]] <- ans[["ee"]] + effects[effNames(c(e, 
                    e), c(l1, l2), nloc)]/(nloc * (nloc - 1)/2)
                }
                else {
                  ans[[paste("e", l1, l2, sep = "")]] <- effects[effNames(c(e, 
                    e), c(l1, l2), nloc)]
                }
            }
        }
    }
    return(ans)
}
