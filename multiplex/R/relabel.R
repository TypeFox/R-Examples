relabel <-
function (S, labels = NULL) 
{
    if (isTRUE(attr(S, "class")[1] == "Semigroup") == FALSE) {
        stop("\"S\" should be an object of a \"Semigroup\" class.")
    }
    ifelse(isTRUE(is.null(labels) == TRUE) == TRUE, labels <- S$st, 
        NA)
    lb <- S$st
    S$st <- labels
    x <- S$S
    ifelse(isTRUE(attr(S, "class")[2] == "symbolic") == TRUE, 
        rownames(x) <- colnames(x) <- labels, NA)
    for (i in which(labels != lb)) {
        x <- replace(x, S$S == as.character(lb[i]), as.character(labels[i]))
    }
    rm(i)
    S$S <- as.data.frame(x)
    rm(x, lb)
    return(S)
}
