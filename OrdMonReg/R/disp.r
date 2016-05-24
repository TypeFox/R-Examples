disp <- function(r, digits){
    res <- format(round(r, digits), nsmall = digits, scientific = FALSE)
    return(res)
    }
