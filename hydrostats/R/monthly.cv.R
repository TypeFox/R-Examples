monthly.cv <- function(flow.ts) {
    month.means <- aggregate(flow.ts[, "Q"], by = list(month = strftime(flow.ts[["Date"]], format = "%m")), mean, na.rm = T)
    month.sd <- sd(month.means[, "x"], na.rm = T)
    
    monthly.cv <- (month.sd/mean(month.means[, "x"])) * 100
    
    return(data.frame(monthly.cv = monthly.cv))
} 
