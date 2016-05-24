convert <-
function (x, labels = NULL, SemigroupClass = FALSE) 
{
    if (isTRUE(attr(x, "class")[1] == "Semigroup") == FALSE) 
        stop("\"x\" should be an object of a \"Semigroup\" class.")
    S <- x$S
    if (isTRUE(attr(x, "class")[2] == "numerical") == TRUE) {
        ifelse(is.null(labels) == FALSE, lbs <- as.vector(labels), 
            lbs <- LETTERS[1:x$ord])
    }
    else if (isTRUE(attr(x, "class")[2] == "symbolic") == TRUE) {
        if (is.null(labels) == FALSE) 
            warning("Semigroup is in 'symbolic' form and 'labels' are therefore ignored.")
        ifelse(is.null(x$st) == TRUE, lbs <- LETTERS[1:x$ord], 
            lbs <- x$st)
    }
    s <- vector()
    for (i in 1:length(as.matrix(S))) {
        if (isTRUE(attr(x, "class")[2] == "symbolic") == TRUE) {
            s[i] <- which(lbs == as.matrix(S)[i])
        }
        else if (isTRUE(attr(x, "class")[2] == "numerical") == 
            TRUE) {
            s[i] <- lbs[which(dimnames(S)[[1]] == as.matrix(S)[i])]
        }
    }
    s <- matrix(s, nrow = nrow(S), ncol = ncol(S))
    if (isTRUE(attr(x, "class")[2] == "numerical") == TRUE) {
        dimnames(s)[[1]] <- dimnames(s)[[2]] <- as.list(lbs)
    }
    else if (isTRUE(attr(x, "class")[2] == "symbolic") == TRUE) {
        dimnames(s)[[1]] <- dimnames(s)[[2]] <- as.list(1:x$ord)
    }
    if (SemigroupClass) {
        lst <- list(ord = nrow(s), st = dimnames(s)[[1]], S = as.data.frame(s))
        ifelse(isTRUE(attr(x, "class")[2] == "numerical") == 
            TRUE, class(lst) <- c("Semigroup", "symbolic"), class(lst) <- c("Semigroup", 
            "numerical"))
        return(lst)
    }
    else {
        return(as.data.frame(s))
    }
}
