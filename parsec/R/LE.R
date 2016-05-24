LE <- function(
    Lambda
) {
    Lambda <- incidence2cover(Lambda)
    n <- nrow(Lambda)
    varnames <- rownames(Lambda)
    if (is.null(varnames)) {
        varnames <- LETTERS[1:n]
        rownames(Lambda)
    }
    colnames(Lambda) <- varnames
    
    lle <- function(cover, mins=NULL) {
        if (dim(cover)[1] < 2)
            return(append(mins, rownames(cover))) 
        minimal <- which(colSums(cover)==0)
        lapply(minimal, function(x) {
            newC <- as.matrix(cover[-x, -x])
            rownames(newC) <- colnames(newC) <- rownames(cover)[-x]
            newMins <- append(mins, rownames(cover)[x])
            lle(newC, newMins)
        })
    }
    
    varords <- as.list(as.data.frame(matrix(unlist(lle(Lambda)), n)))
    varords <- lapply(varords, as.character)
    names(varords) <- paste("LE", 1:length(varords), sep="")
    return(varords)
}