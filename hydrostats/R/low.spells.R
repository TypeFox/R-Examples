low.spells <- function(flow.ts, quant = 0.1, threshold = NULL, duration = T, volume = T, plot = T, ann.stats = T, ann.stats.only = F, hydro.year = FALSE) {
    gauge <- deparse(substitute(flow.ts))
    
    if (hydro.year == TRUE) {
        print("Returning results based on hydrologic year")
        flow.ts <- hydro.year(flow.ts, year = "hydro")
        record.year <- flow.ts[["hydro.year"]]
    } else {
        record.year <- strftime(flow.ts[["Date"]], format = "%Y")
        flow.ts <- data.frame(flow.ts, hydro.year = record.year)
    }
    
    n.years <- nlevels(as.factor(record.year))
    
    if (ann.stats == T) {
        # calculate annual minimum flow
        flow.ts.comp <- na.omit(flow.ts)
        
        n.days <- tapply(flow.ts.comp[["Q"]], flow.ts.comp[["hydro.year"]], length)
        n.most.days <- which(n.days >= 350)
        if (length(n.most.days) == 0) {
            ann.mins.mean <- NA
            cv.min.ann <- NA
            avg.min.day <- data.frame(mean.doy = NA, sd.doy = NA)
            min.min.ann.duration <- NA
            max.min.ann.duration <- NA
            avg.min.ann.duration <- NA
            
            
        } else {
            
            
            flow.ts.comp <- flow.ts.comp[which(flow.ts.comp[["hydro.year"]] %in% names(n.most.days)), ]
            flow.ts.comp[["hydro.year"]] <- factor(flow.ts.comp[["hydro.year"]])
            record.year <- flow.ts.comp[["hydro.year"]]
            n.years <- nlevels(as.factor(record.year))
            
            
            ann.mins <- aggregate(flow.ts.comp["Q"], flow.ts.comp["hydro.year"], min, na.rm = T)
            
            ann.min.days <- merge(ann.mins, flow.ts.comp)
            
            
            ann.mins.mean <- mean(ann.mins[["Q"]], na.rm = T)
            ann.mins.sd <- sd(ann.mins[["Q"]], na.rm = T)
            
            
            avg.ann.min.days <- aggregate(ann.min.days["Date"], ann.min.days["hydro.year"], function(x) t(day.dist(x)))
            
            avg.min.day <- day.dist(days = avg.ann.min.days[["Date"]][[1]], years = avg.ann.min.days[["hydro.year"]])
            
            min.ann.flow.threshold <- max(ann.mins$Q, na.rm = T)
            
            ann.min.spells <- ifelse(flow.ts.comp[["Q"]] <= min.ann.flow.threshold, 1, 0)
            # ann.max.spell.runs <- rle(ann.max.spells) now done seperately for each year below.
            
            ann.runs <- tapply(ann.min.spells, flow.ts.comp[["hydro.year"]], rle)
            ann.runs.cum <- lapply(ann.runs, function(x) cumsum(x$lengths))
            
            
            ann.min.index.all.events <- tapply(flow.ts.comp[["Q"]], flow.ts.comp[["hydro.year"]], function(x) which(x == min(x)))
            
            if (max(sapply(ann.min.index.all.events, length)) > 1) {
                ann.min.event.rles <- lapply(ann.min.index.all.events, function(x) rle(diff(x)))
                idx <- !(sapply(ann.min.event.rles, function(x) length(x$lengths)))
                ann.min.event.rles[idx] <- lapply(ann.min.event.rles[idx], function(x) {
                  x$lengths <- 0
                  x$values <- 1
                  x
                })
                
                ann.min.events.max.dur <- unlist(lapply(ann.min.event.rles, function(x) max(x$lengths[which(x$values == 1)]) + 1))  #add 1 because diff function returns index-1
                ann.min.events.max.dur <- ann.min.events.max.dur[is.finite(ann.min.events.max.dur)]
            } else {
                ann.min.rle.event.indices <- sapply(ann.runs, function(x) which(x$values == 1))
                
                ann.min.rle.event.indices <- mapply(function(x, y) x[y], ann.runs.cum, ann.min.rle.event.indices)
                ann.min.rle.run.indices <- lapply(ann.runs, function(x) which(x$values == 1))
                
                smallest.event.index <- mapply(function(x, y) which.min(abs(x - y)), ann.min.index.all.events, ann.min.rle.event.indices)
                final.event.index <- mapply(function(x, y) x[y], ann.min.rle.run.indices, smallest.event.index)
                
                ann.min.events.max.dur <- unlist(mapply(function(x, y) x$lengths[y], ann.runs, final.event.index, SIMPLIFY = TRUE))
                
            }
            
            
            min.min.ann.duration <- min(ann.min.events.max.dur, na.rm = T)
            avg.min.ann.duration <- mean(ann.min.events.max.dur, na.rm = T)
            max.min.ann.duration <- max(ann.min.events.max.dur, na.rm = T)
            
            cv.min.ann <- (ann.mins.sd/ann.mins.mean) * 100
        }
    }
    
    if (ann.stats.only == T) {
        
        
        
        return(data.frame(avg.min.ann = ann.mins.mean, cv.min.ann = cv.min.ann, ann.min.timing = avg.min.day[, "mean.doy"], ann.min.timing.sd = avg.min.day[, "sd.doy"], ann.min.min.dur = min.min.ann.duration, 
            ann.min.avg.dur = avg.min.ann.duration, ann.min.max.dur = max.min.ann.duration))
        
    } else {
        
        if (!is.null(threshold)) {
            
            
            flow.threshold <- threshold
            
        } else {
            
            
            # average spell characteristics
            flow.threshold <- quantile(flow.ts[, "Q"], quant, na.rm = T)
            names(flow.threshold) <- NULL
            
        }
        low.flows <- ifelse(flow.ts[, "Q"] <= flow.threshold, 1, 0)
        low.flow.av <- mean(flow.ts[which(low.flows == 1), "Q"])
        low.flow.sd <- sd(flow.ts[which(low.flows == 1), "Q"])
        
        
        low.flow.runs <- rle(low.flows)
        
        
        low.spell.days <- as.numeric(strftime(flow.ts[which(low.flows == 1), "Date"], format = "%j"))
        
        good.low.flow.runs <- which(!is.na(low.flow.runs$values))
        flow.runs.values <- low.flow.runs$values[good.low.flow.runs]
        flow.runs.lengths <- low.flow.runs$lengths[good.low.flow.runs]
        
        low.spell.frequency <- length(flow.runs.values[flow.runs.values == 1])/n.years
        
        
        
        
        
        
        if (duration == TRUE) {
            min.duration <- min(low.flow.runs$lengths[which(low.flow.runs$values == 1)], na.rm = T)
            avg.duration <- mean(low.flow.runs$lengths[which(low.flow.runs$values == 1)], na.rm = T)
            med.duration <- median(low.flow.runs$lengths[which(low.flow.runs$values == 1)], na.rm = T)
            max.duration <- max(low.flow.runs$lengths[which(low.flow.runs$values == 1)], na.rm = T)
            sd.duration <- sd(low.flow.runs$lengths[which(low.flow.runs$values == 1)], na.rm = T)
            cv.duration <- sd.duration/avg.duration
        }
        
        if (volume == TRUE) {
            spell.factor <- rep(seq_along(low.flow.runs$lengths), times = low.flow.runs$lengths)
            spells <- split(flow.ts[, "Q"], spell.factor)
            spell.volumes <- flow.ts[, "Q"]
            spell.volumes <- sapply(spells, sum)
            spell.volumes <- spell.volumes[which(low.flow.runs$values == 1)]
            
        }
        
        
        if (plot == TRUE) {
            plot(flow.ts[, "Date"], flow.ts[, "Q"], type = "l", main = gauge, xlab = "Date", ylab = "Q")
            
            points(flow.ts[which(low.flows == 1), "Date"], flow.ts[which(low.flows == 1), "Q"], col = "red", cex = 0.25)
            
            abline(h = flow.threshold)
        }
        
        
    }
    
    if (ann.stats == F) {
        return(data.frame(low.spell.threshold = flow.threshold, min.low.spell.duration = min.duration, avg.low.spell.duration = avg.duration, med.low.spell.duration = med.duration, max.low.duration = max.duration, 
            low.spell.freq = low.spell.frequency))
        
    } else {
        
        return(data.frame(low.spell.threshold = flow.threshold, min.low.spell.duration = min.duration, avg.low.spell.duration = avg.duration, med.low.spell.duration = med.duration, max.low.duration = max.duration, 
            low.spell.freq = low.spell.frequency, avg.min.ann = ann.mins.mean, cv.min.ann = cv.min.ann, ann.min.timing = avg.min.day[, "mean.doy"], ann.min.timing.sd = avg.min.day[, "sd.doy"], 
            ann.min.min.dur = min.min.ann.duration, ann.min.avg.dur = avg.min.ann.duration, ann.min.max.dur = max.min.ann.duration))
    }
    
    
} 
