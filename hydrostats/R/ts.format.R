ts.format <- function(x, format = "%d/%m/%Y", cols = c(1, 2)) {
    names(x)[cols] <- c("Date", "Q")
    x[["Date"]] <- strptime(x[["Date"]], format = format)
    x[["Date"]] <- as.POSIXct(x[["Date"]])
    return(x)
} 
