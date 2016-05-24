Colwells <- function(flow.ts, fn = "mean", boundaries = "transform", s = 11, base = 2, from = 0.5, by = 0.25, indices.only = FALSE) {
    fn <- match.fun(fn)
    
    
    # fn<-summary.stat
    
    flow.ts$month <- factor(strftime(flow.ts[, "Date"], format = "%m"))
    flow.ts$year <- factor(strftime(flow.ts[, "Date"], format = "%Y"))
    
    flow.ts.monthly <- aggregate(Q ~ month + year, flow.ts, fn, na.rm = TRUE)
    
    
    
    if (boundaries == "transform") {
        
        flow.ts.monthly$Q <- log10(flow.ts.monthly$Q + 1)
        flow.ts.monthly$flow.class <- cut(flow.ts.monthly$Q, s, right = FALSE, include.lowest = TRUE)
        flow.table <- with(flow.ts.monthly, table(flow.class, month))
        pbreaks <- "see Table"
        
        
    } else {
        
        if (boundaries == "log_class_size") {
            
            breaks <- base^(seq(1:s) - 2)
            breaks <- c(-Inf, breaks, Inf)
            pbreaks <- breaks
            flow.ts.monthly$flow.class <- cut(flow.ts.monthly$Q, breaks)
            flow.table <- with(flow.ts.monthly, table(flow.class, month))
            
        } else {
            
            if (boundaries == "weighted_log_class_size") {
                low_exp <- ceiling(0 - (s - 1)/2)
                breaks <- base^seq(low_exp, length.out = s - 1)
                breaks <- c(-Inf, breaks, Inf)
                pbreaks <- breaks
                breaks <- fn(flow.ts.monthly$Q) * breaks
                
                
                flow.ts.monthly$flow.class <- cut(flow.ts.monthly$Q, breaks, right = FALSE, include.lowest = TRUE)
                flow.table <- with(flow.ts.monthly, table(flow.class, month))
                
                
            } else {
                
                if (boundaries == "equal") {
                  flow.ts.monthly$flow.class <- cut(flow.ts.monthly$Q, s, right = FALSE, include.lowest = TRUE)
                  flow.table <- with(flow.ts.monthly, table(flow.class, month))
                  pbreaks <- "see Table"
                } else {
                  
                  if (boundaries == "Gan") {
                    Q <- fn(flow.ts.monthly$Q, na.rm = TRUE)
                    breaks <- seq(from = from, by = by, length.out = s - 1)
                    pbreaks <- breaks
                    breaks <- Q * breaks
                    breaks <- c(-Inf, breaks, Inf)
                    flow.ts.monthly$flow.class <- cut(flow.ts.monthly$Q, breaks, right = FALSE, include.lowest = TRUE)
                    flow.table <- with(flow.ts.monthly, table(flow.class, month))
                  }
                }
            }
        }
    }
    # X<-margin.table(flow.table, 2)
    X <- apply(flow.table, 2, sum, na.rm = T)
    # Y<-margin.table(flow.table, 1)
    Y <- apply(flow.table, 1, sum, na.rm = T)
    Z <- sum(flow.table, na.rm = TRUE)
    
    HX <- -1 * sum((X/Z) * log(X/Z), na.rm = TRUE)
    HY <- -1 * sum((Y/Z) * log(Y/Z), na.rm = TRUE)
    HXY <- -1 * sum((flow.table/Z) * log(flow.table/Z), na.rm = TRUE)
    
    P <- round(1 - (HXY - HX)/log(s), 2)
    C <- round(1 - HY/log(s), 2)
    M <- round((HX + HY - HXY)/log(s), 2)
    
    
    if (indices.only == TRUE) {
        
        return(data.frame(P = P, C = C, M = M, CP = C/P, MP = M/P))
        
    } else {
        
        return(list(breaks = pbreaks, flow.table = flow.table, P = P, C = C, M = M, CP = C/P, MP = M/P))
    }
} 
