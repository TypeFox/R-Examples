dichot <-
function (x, c = 1, diag.incl = TRUE) 
{
    if (isTRUE(diag.incl == TRUE) == FALSE) {
        if (is.na(dim(x)[3]) == TRUE | isTRUE(dim(x)[3] == 1) == 
            TRUE) {
            dg <- diag(x)
        }
        else {
            dg <- list()
            for (i in 1:dim(x)[3]) dg[[i]] <- diag(x[, , i])
        }
    }
    else {
        dg <- NULL
    }
    if (isTRUE(0 >= c) == TRUE) {
        x <- replace(x, x >= c, 1L)
        x <- replace(x, x < c, 0L)
    }
    else {
        x <- replace(x, x < c, 0L)
        x <- replace(x, x >= c, 1L)
    }
    if (isTRUE(is.null(dg)) == FALSE) {
        if (is.na(dim(x)[3]) == TRUE | isTRUE(dim(x)[3] == 1) == 
            TRUE) {
            diag(x) <- dg
        }
        else {
            for (i in 1:dim(x)[3]) diag(x[, , i]) <- dg[[i]]
        }
        x
    }
    else {
        x
    }
}
