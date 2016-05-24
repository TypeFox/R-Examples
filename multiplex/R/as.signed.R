as.signed <-
function (x, labels = NULL) 
{
    if (is.array(x) == FALSE) 
        stop("Data must be an array")
    if (is.na(dim(x)[3]) == FALSE) {
        sm <- x[, , 1]
        warning("Take the 1st dim. in 'x' only.")
    }
    else {
        sm <- x
    }
    if (is.null(labels) == FALSE) {
        ifelse(isTRUE(length(labels) == dim(sm)[1]) == TRUE, 
            NA, labels <- 1:dim(sm)[1])
        rownames(sm) <- colnames(sm) <- labels
    }
    else if (is.null(dimnames(sm)[1]) == TRUE) {
        rownames(sm) <- colnames(sm) <- 1:dim(sm)[1]
    }
    else {
        NA
    }
    val <- levels(factor(sm))
    lst <- list(val = noquote(levels(stats::reorder(val, length(val):1))), 
        s = sm)
    class(lst) <- "Signed"
    return(lst)
}
