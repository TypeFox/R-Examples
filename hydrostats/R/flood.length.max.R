flood.length.max <- function(flow.ts, threshold, ind.days = 5) {
    
    high.flows <- ifelse(flow.ts[["Q"]] > threshold, 1, 0)
    
    high.flow.runs <- rle(high.flows)
    too.short <- which(high.flow.runs$lengths < ind.days & high.flow.runs$values == 0)
    spell.factor <- rep(seq_along(high.flow.runs$lengths), times = high.flow.runs$lengths)
    add.to.spell <- which(spell.factor %in% too.short)
    high.flows[add.to.spell] <- 1
    
    
    high.flow.runs <- rle(high.flows)
    if (length(which(high.flow.runs$values == 1)) == 0) {
        
        max.duration <- 0
    } else {
        max.duration <- max(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
    }
    
    max.duration <- data.frame(max.duration = max.duration)
    
    return(max.duration)
} 
