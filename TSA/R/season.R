`season` <-
function (x, labels) 
{
    if (!is.ts(x) || frequency(x) == 1) 
        stop("x has to be a time series with freq >1")
    if (missing(labels)) {
        default = paste("Season-", 1:frequency(x), sep = "")
        labels = switch(as.character(frequency(x)), "4" = c("1Q", 
            "2Q", "3Q", "4Q"), "7" = c("Monday", "Tuesday", "Wednesday", 
            "Thursday", "Friday", "Saturday", "Sunday"), "12" = c("January", 
            "February", "March", "April", "May", "June", "July", 
            "August", "September", "October", "November", "December"), 
            default)
    }
    t1 = start(x)[2]
    if (t1 > 1) 
        y = rep(c(t1:frequency(x), 1:(t1 - 1)), ceiling(length(x)/frequency(x)))[1:length(x)]
    else y = rep((1:frequency(x)), ceiling(length(x)/frequency(x)))
    y = factor(y[seq(x)], labels = labels)
    invisible(y)
}

