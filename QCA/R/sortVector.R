 ## function to sort a vector of strings according to their length

`sortVector` <- 
function(x, collapse = "") {
    
    strx <- strsplit(x, split=ifelse(collapse == "", "", "\\*"))
    strings <- NULL
    lengths <- unlist(lapply(strx, length))
    unique.lengths <- sort(unique(lengths))
    for (i in seq(length(unique.lengths))) {
        strings <- c(strings, sort(x[which(lengths == unique.lengths[i])]))
    }
    return(strings)
}

