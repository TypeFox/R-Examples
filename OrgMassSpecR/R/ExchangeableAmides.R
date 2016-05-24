ExchangeableAmides <- function(sequence) {
    n <- length(sequence)
    x <- vector(mode = "numeric", length = n)
    for(i in 1:n) {
        seq_vector <- strsplit(as.character(sequence[i]), split = "")[[1]]
        x[i] <- length(na.omit(sub("P", NA, seq_vector))) - 1
    }
    return(x)	
}
