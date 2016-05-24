startingValuesNothing <-
function (nloc, max.level = 2, max.dom = 2, e.unique = FALSE) 
{
    ans <- list()
    for (l1 in 1:nloc) {
        ans[[paste("a", l1, sep = "")]] <- 0
        if (max.dom > 0) {
            ans[[paste("d", l1, sep = "")]] <- 0
        }
    }
    if ((max.level > 1) && (nloc > 1)) {
        if (e.unique) {
            ans[["ee"]] <- 0
        }
        else {
            for (l1 in 1:(nloc - 1)) {
                for (l2 in (l1 + 1):nloc) {
                  ans[[paste("e", l1, l2, sep = "")]] <- 0
                }
            }
        }
    }
    return(ans)
}
