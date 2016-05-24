dhc <-
function (xd, prsep = ", ") 
{
    if (isTRUE(length(xd) != 0L) == TRUE && isTRUE(is.na(xd)) == 
        FALSE) {
        Ltd <- FALSE
        ifelse(isTRUE(is.list(xd)) == TRUE, Ltd <- TRUE, xd <- as.list(xd))
        Xd <- list()
        length(Xd) <- length(xd)
        for (i in 1:length(xd)) {
            if (isTRUE(length(xd[[i]]) != 0L) == TRUE) {
                tmpd <- as.list(xd[[i]])
                for (j in 1:length(xd[[i]])) Xd[[i]] <- append(Xd[[i]], 
                  strsplit(tmpd[[j]], prsep)[[1]])
            }
        }
        attr(Xd, "names") <- attr(xd, "names")
        ifelse(isTRUE(Ltd) == TRUE, return(Xd), return(unlist(Xd)))
    }
    else {
        xd
    }
}
