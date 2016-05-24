CTF <- function(flow.ts, threshold = 0.1) {
    
    flow.ts <- na.omit(flow.ts)
    
    if (length(flow.ts["Q"]) == 0 | sum(flow.ts[, "Q"]) == 0) {
        
        return(data.frame(p.CTF = 1, avg.CTF = NA, med.CTF = NA, min.CTF = NA, max.CTF = NA))
    }
    
    
    if (sum(flow.ts["Q"] <= threshold) < 1) {
        
        return(data.frame(p.CTF = 0, avg.CTF = 0, med.CTF = 0, min.CTF = 0, max.CTF = 0))
    } else {
        
        CTF.ts <- ifelse(flow.ts["Q"] > threshold, 1, 0)
        
        CTF.runs <- rle(CTF.ts[, "Q"])
        p.CTF <- sum(CTF.ts == 0)/length(CTF.ts)
        avg.CTF <- mean(CTF.runs$lengths[which(CTF.runs$values == 0)], na.rm = T)
        med.CTF <- median(CTF.runs$lengths[which(CTF.runs$values == 0)], na.rm = T)
        min.CTF <- min(CTF.runs$lengths[which(CTF.runs$values == 0)], na.rm = T)
        max.CTF <- max(CTF.runs$lengths[which(CTF.runs$values == 0)], na.rm = T)
        
        
    }
    return(data.frame(p.CTF = p.CTF, avg.CTF = avg.CTF, med.CTF = med.CTF, min.CTF = min.CTF, max.CTF = max.CTF))
    
    
} 
