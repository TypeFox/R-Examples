ann.cv <- function(flow.ts) {
    flow.ts$year <- strftime(flow.ts$Date, format = "%Y")
    ann.stats <- aggregate(flow.ts["Q"], flow.ts["year"], mean, na.rm = T)
    
    
    out <- sd(ann.stats$Q, na.rm = T)/mean(ann.stats$Q, na.rm = T) * 100
    
    return(data.frame(ann.cv = out))
} 
