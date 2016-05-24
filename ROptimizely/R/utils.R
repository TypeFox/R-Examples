empty_list_check <- function(x) {
    if (length(x) == 0 & is.list(x)) {
        return(NA)
    } else if (length(x) > 1) {
        return(paste(x, collapse = ","))
    } else {
        return(x)
    }
} 
