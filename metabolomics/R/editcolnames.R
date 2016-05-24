editcolnames <- function(y)
{
    colnames(y) <- if (
        length(
            grep("^X[\\d]", colnames(y), perl=TRUE)
        ) != 0
    ) {
        gsub("^X([\\d].*)", "\\1", colnames(y), perl=TRUE)
    } else {
        colnames(y)
    }
    
    return(y)
}
