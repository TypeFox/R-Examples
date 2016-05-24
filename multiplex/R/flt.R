flt <-
function (e, PO, rclos = TRUE) 
{
    if (isTRUE(e > nrow(PO)) == TRUE) 
        stop("'e' is greater than the size of the partial order")
    pfl <- vector()
    for (i in 1:nrow(PO)) {
        if (isTRUE(PO[e, i] == 1) == TRUE && isTRUE(PO[i, e] == 
            0) == TRUE) {
            pfl <- append(pfl, i)
        }
        else {
            NA
        }
    }
    rm(i)
    if (rclos) {
        pfl <- append(pfl, e)
    }
    return(pfl)
}
