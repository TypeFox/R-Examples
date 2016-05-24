partial.series <- function(flow.ts, ari = 2, ind.days = 7, duration = T, plot = F, volume = T, series = FALSE) {
    gauge <- deparse(substitute(flow.ts))
    
    record.year <- strftime(flow.ts[["Date"]], format = "%Y")
    flow.ts <- data.frame(flow.ts, hydro.year = record.year)
    flow.ts.comp <- na.omit(flow.ts)
    n.days <- tapply(flow.ts.comp[["Q"]], flow.ts.comp[["hydro.year"]], length)
    n.most.days <- which(n.days >= 350)
    flow.ts.comp <- flow.ts.comp[which(flow.ts.comp[["hydro.year"]] %in% names(n.most.days)), ]
    flow.ts.comp[["hydro.year"]] <- factor(flow.ts.comp[["hydro.year"]])
    record.year <- flow.ts.comp[["hydro.year"]]
    
    n.years <- length(n.most.days)
    
    
    if (ari > n.years) {
        print("warning(time-series is shorter than ari. Please provide a smaller ari")
        if (volume == FALSE) {
            return(data.frame(ari = ari, n.years = n.years, n.events = NA, flow.threshold = NA, avg.duration = NA, max.duration = NA))
        } else {
            return(data.frame(ari = ari, n.years = n.years, n.events = NA, flow.threshold = NA, avg.duration = NA, max.duration = NA, med.spell.volume = NA))
        }
        
        
    }
    n.events <- data.frame(n.events = ceiling(n.years/ari))
    p.series <- vector("list", length = n.events$n.events)
    
    rising <- data.frame(rising = flow.ts[2:nrow(flow.ts), "Q"] - flow.ts[1:nrow(flow.ts) - 1, "Q"])
    falling <- data.frame(falling = flow.ts[3:nrow(flow.ts), "Q"] - flow.ts[2:nrow(flow.ts) - 2, "Q"])
    
    peak.search <- data.frame(flow.ts, rising = c(NA, rising[["rising"]]), falling = c(falling[["falling"]], NA, NA))
    peaks <- flow.ts[which(peak.search[["rising"]] > 0 & peak.search[["falling"]] < 0), ]
    
    peaks.ord <- peaks[order(peaks[["Q"]], decreasing = T), ]
    
    p.series[[1]] <- data.frame(peaks.ord[1, ])
    
    i = 1
    
    while (i < n.events$n.events) {
        i <- i + 1
        dif.time.test <- difftime(peaks.ord$Date[1], peaks.ord$Date)
        
        peaks.ord <- peaks.ord[which(abs(dif.time.test) > (ind.days * 24 * 60 * 60)), ]
        
        p.series[[i]] <- data.frame(peaks.ord[1, ])
        
        if (is.na(peaks.ord[1, "Q"]) == T) 
            NA
    }
    p.series <- do.call("rbind", p.series)
    
    p.series$event.rank <- seq(1:nrow(p.series))
    
    flow.threshold <- tail(p.series[["Q"]], 1)
    
    if (plot == TRUE) {
        plot(flow.ts[["Date"]], flow.ts[["Q"]], type = "l", main = gauge, xlab = "Date", ylab = "Q")
        
        points(p.series$Date, p.series$Q, col = "red", cex = 0.25)
        abline(h = (tail(p.series["Q"], 1) - 1))
    }
    
    high.flows <- ifelse(flow.ts[["Q"]] >= flow.threshold, 1, 0)
    high.flow.runs <- rle(high.flows)
    
    if (duration == TRUE) {
        avg.duration <- mean(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
        max.duration <- max(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
    }
    if (volume == TRUE) {
        spell.factor <- rep(seq_along(high.flow.runs$lengths), times = high.flow.runs$lengths)
        spells <- split(flow.ts[["Q"]], spell.factor)
        spell.volumes <- flow.ts[["Q"]]
        spell.volumes <- sapply(spells, sum)
        spell.volumes.below.threshold <- sapply(spells, length) * flow.threshold
        spell.volumes <- spell.volumes[which(high.flow.runs$values == 1)] - spell.volumes.below.threshold[which(high.flow.runs$values == 1)]
        
        
        if (series == TRUE) {
            return(list(p.series = p.series, ari = ari, n.years = n.years, n.events = n.events, flow.threshold = flow.threshold, avg.duration = avg.duration, max.duration = max.duration, med.spell.volume = median(spell.volumes)))
        } else {
            return(data.frame(ari = ari, n.years = n.years, n.events = n.events, flow.threshold = flow.threshold, avg.duration = avg.duration, max.duration = max.duration, med.spell.volume = median(spell.volumes)))
        }
    } else {
        if (series == TRUE) {
            return(list(p.series = p.series, ari = ari, n.years = n.years, n.events = n.events, flow.threshold = flow.threshold, avg.duration = avg.duration, max.duration = max.duration))
        } else {
            return(data.frame(ari = ari, n.years = n.years, n.events = n.events, flow.threshold = flow.threshold, avg.duration = avg.duration, max.duration = max.duration))
        }
    }
} 
