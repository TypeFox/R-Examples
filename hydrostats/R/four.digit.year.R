four.digit.year <- function(x, year = 1968) {
    n <- as.numeric(strftime(x, format = "%y"))%%100
    Y <- ifelse(n > year%%100, 1900 + n, 2000 + n)
    return(Y)
} 
