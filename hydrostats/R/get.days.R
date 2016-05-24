get.days <- function(year) {
    days <- c()
    d1 <- as.Date(paste(year, "-01-01", sep = ""))
    d2 <- as.Date(paste(year, "-12-31", sep = ""))
    day <- difftime(d2, d1) + 1
    days <- rbind(day, days)
    return(as.vector(days))
} 
