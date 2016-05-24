daily.cv <- function(flow.ts) {
    
    daily.cv <- (sd(flow.ts$Q, na.rm = T)/mean(flow.ts$Q, na.rm = T)) * 100
    
    return(data.frame(daily.cv = daily.cv))
} 
