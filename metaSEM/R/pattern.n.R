pattern.na <- function(x, show.na=TRUE) {
    out <- Reduce("+", lapply(x, is.na))
    if (show.na) out else length(x)-out
} 
