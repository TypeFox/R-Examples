high.spells <- function(flow.ts, quant = 0.9, threshold = NULL, ind.days = 5, duration = TRUE, volume = TRUE, plot = TRUE, ignore.zeros = FALSE, ctf.threshold = 0.1, ann.stats = TRUE, ann.stats.only = FALSE, 
    inter.flood = FALSE, hydro.year = FALSE) {
    
    gauge <- deparse(substitute(flow.ts))
    
    
    if (hydro.year == TRUE) {
        print("Returning results based on hydrologic year")
        flow.ts <- hydro.year(flow.ts, hydro.year = "hydro")
        record.year <- flow.ts[, "hydro.year"]
    } else {
        record.year <- strftime(flow.ts[["Date"]], format = "%Y")
        flow.ts <- data.frame(flow.ts, hydro.year = record.year)
    }
    n.years <- nlevels(as.factor(record.year))
    
    if (ann.stats == T) {
        
        flow.ts.comp <- na.omit(flow.ts)
        
        n.days <- tapply(flow.ts.comp[["Q"]], flow.ts.comp[["hydro.year"]], length)
        n.most.days <- which(n.days >= 350)
        
        if (length(n.most.days) == 0) {
            ann.maxs.mean <- NA
            ann.maxs.sd <- NA
            avg.max.day <- data.frame(mean.doy = NA, sd.doy = NA)
            min.max.ann.duration <- NA
            avg.max.ann.duration <- NA
            max.max.ann.duration <- NA
            flood.skewness <- NA
            cv.ann.duration <- NA
            cv.max.ann <- NA
        } else {
            flow.ts.comp <- flow.ts.comp[which(flow.ts.comp[["hydro.year"]] %in% names(n.most.days)), ]
            flow.ts.comp[["hydro.year"]] <- factor(flow.ts.comp[["hydro.year"]])
            record.year <- flow.ts.comp[["hydro.year"]]
            n.years <- nlevels(as.factor(record.year))
            
            
            ann.maxs <- aggregate(flow.ts.comp["Q"], flow.ts.comp["hydro.year"], max, na.rm = T)
            
            mean.ann <- data.frame(mean = mean(flow.ts.comp[["Q"]], na.rm = T))
            
            
            ann.max.days <- merge(ann.maxs, flow.ts.comp)
            
            
            ann.maxs.mean <- mean(ann.maxs[["Q"]], na.rm = T)
            
            ann.maxs.sd <- sd(ann.maxs[["Q"]], na.rm = T)
            
            avg.ann.max.days <- aggregate(ann.max.days["Date"], ann.max.days["hydro.year"], function(x) t(day.dist(x)))
            
            avg.max.day <- day.dist(days = avg.ann.max.days[["Date"]][[1]], years = avg.ann.max.days[["hydro.year"]])
            
            
            
            max.ann.flow.threshold <- min(ann.maxs$Q, na.rm = T)
            
            ann.max.spells <- ifelse(flow.ts.comp[["Q"]] >= max.ann.flow.threshold, 1, 0)
            # ann.max.spell.runs <- rle(ann.max.spells) now done seperately for each year below.
            
            ann.runs <- tapply(ann.max.spells, flow.ts.comp[["hydro.year"]], rle)
            ann.runs.cum <- lapply(ann.runs, function(x) cumsum(x$lengths))
            
            ann.max.index.all.events <- tapply(flow.ts.comp[["Q"]], flow.ts.comp[["hydro.year"]], function(x) which(x == max(x)))
            
            if (max(sapply(ann.max.index.all.events, length)) > 1) {
                ann.max.event.rles <- lapply(ann.max.index.all.events, function(x) rle(diff(x)))
                idx <- !(sapply(ann.max.event.rles, function(x) length(x$lengths)))
                ann.max.event.rles[idx] <- lapply(ann.max.event.rles[idx], function(x) {
                  x$lengths <- 0
                  x$values <- 1
                  x
                })
                
                ann.max.events.max.dur <- unlist(lapply(ann.max.event.rles, function(x) max(x$lengths[which(x$values == 1)], na.rm=T) + 1))  #add 1 because diff function returns index-1
                ann.max.events.max.dur <- ann.max.events.max.dur[is.finite(ann.max.events.max.dur)]
                
            } else {
                
                ann.max.rle.event.indices <- sapply(ann.runs, function(x) which(x$values == 1))  #list of vectors
                
                ann.max.rle.event.indices <- mapply(function(x, y) x[y], ann.runs.cum, ann.max.rle.event.indices)  #list
                
                ann.max.rle.run.indices <- lapply(ann.runs, function(x) which(x$values == 1))  #list
                
                largest.event.index <- mapply(function(x, y) which.min(abs(x - y)), ann.max.index.all.events, ann.max.rle.event.indices)  #vector
                
                final.event.index <- mapply(function(x, y) x[y], ann.max.rle.run.indices, largest.event.index)  #vector
                
                ann.max.events.max.dur <- unlist(mapply(function(x, y) x$lengths[y], ann.runs, final.event.index, SIMPLIFY = TRUE))  #vector
                
            }
            
            
            
            min.max.ann.duration <- min(ann.max.events.max.dur, na.rm = T)
            avg.max.ann.duration <- mean(ann.max.events.max.dur, na.rm = T)
            max.max.ann.duration <- max(ann.max.events.max.dur, na.rm = T)
            cv.max.ann <- (ann.maxs.sd/ann.maxs.mean) * 100
            flood.skewness <- ann.maxs.mean/mean.ann[["mean"]]
            cv.ann.duration <- (sd(ann.max.events.max.dur))/mean(ann.max.events.max.dur) * 100
            
        }
    }
    if (ann.stats.only == T) {
        
        return(data.frame(avg.max.ann = ann.maxs.mean, cv.max.ann = cv.max.ann, ann.max.timing = avg.max.day[["mean.doy"]], ann.max.timing.sd = avg.max.day[["sd.doy"]], flood.skewness = flood.skewness, 
            ann.max.min.dur = min.max.ann.duration, ann.max.avg.dur = avg.max.ann.duration, ann.max.max.dur = max.max.ann.duration, ann.max.cv.dur = cv.ann.duration))
        
    } else {
        
        if (!is.null(threshold)) {
            
            
            flow.threshold <- threshold
            
        } else {
            
            
            if (ignore.zeros == T) {
                
                flow.threshold <- quantile(flow.ts[which(flow.ts[, "Q"] > ctf.threshold), "Q"], quant, na.rm = T)
                names(flow.threshold) <- NULL
            } else {
                flow.threshold <- quantile(flow.ts[, "Q"], quant, na.rm = T)
                names(flow.threshold) <- NULL
            }
            
        }
        
        rise.fall <- c(NA, flow.ts[2:nrow(flow.ts), "Q"] - flow.ts[1:nrow(flow.ts) - 1, "Q"])
        
        
        
        high.flows <- ifelse(flow.ts[, "Q"] >= flow.threshold, 1, 0)
        
        if (ind.days > 0) {
            
            high.flow.runs <- rle(high.flows)
            too.short <- which(high.flow.runs$lengths < ind.days & high.flow.runs$values == 0)
            spell.factor <- rep(seq_along(high.flow.runs$lengths), times = high.flow.runs$lengths)
            add.to.spell <- which(spell.factor %in% too.short)
            high.flows[add.to.spell] <- 1
        }
        
        
        high.flow.runs <- rle(high.flows)
        
        
        
        high.flow.av <- mean(flow.ts[which(high.flows == 1), "Q"], na.rm = T)
        high.flow.sd <- sd(flow.ts[which(high.flows == 1), "Q"], na.rm = T)
        
        high.flow.rf <- rise.fall[which(high.flows == 1)]
        mean.rise <- abs(mean(high.flow.rf[high.flow.rf > 0], na.rm = T))
        mean.fall <- abs(mean(high.flow.rf[high.flow.rf < 0], na.rm = T))
        
        good.high.flow.runs <- which(!is.na(high.flow.runs$values))
        flow.runs.values <- high.flow.runs$values[good.high.flow.runs]
        flow.runs.lengths <- high.flow.runs$lengths[good.high.flow.runs]
        
        n.events <- length(flow.runs.values[flow.runs.values == 1])
        flood.frequency <- n.events/n.years
        
        
        
        
        if (duration == TRUE) {
            
            min.duration <- min(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
            avg.duration <- mean(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
            med.duration <- quantile(high.flow.runs$lengths[which(high.flow.runs$values == 1)], 0.5, na.rm = T, names = F)
            max.duration <- max(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
            sd.duration <- sd(high.flow.runs$lengths[which(high.flow.runs$values == 1)], na.rm = T)
            cv.duration <- sd.duration/avg.duration * 100
        }
        
        if (inter.flood == TRUE) {
            
            avg.interval <- mean(high.flow.runs$lengths[which(high.flow.runs$values == 0)], na.rm = T)
            min.interval <- min(high.flow.runs$lengths[which(high.flow.runs$values == 0)], na.rm = T)
            max.interval <- max(high.flow.runs$lengths[which(high.flow.runs$values == 0)], na.rm = T)
            
            sd.interval <- sd(high.flow.runs$lengths[which(high.flow.runs$values == 0)], na.rm = T)
            cv.interval <- sd.interval/avg.interval * 100
            return(data.frame(high.spell.threshold = flow.threshold, average.interval = avg.interval, min.interval = min.interval, max.interval = max.interval))
        }
        
        if (volume == TRUE) {
            spell.factor <- rep(seq_along(high.flow.runs$lengths), times = high.flow.runs$lengths)
            spells <- split(flow.ts[, "Q"], spell.factor)
            spell.volumes <- flow.ts[, "Q"]
            spell.volumes <- sapply(spells, sum)
            spell.volumes.below.threshold <- sapply(spells, length) * flow.threshold
            spell.volumes <- spell.volumes[which(high.flow.runs$values == 1)] - spell.volumes.below.threshold[which(high.flow.runs$values == 1)]
            
        }
        
        
        if (plot == TRUE) {
            plot(flow.ts[["Date"]], flow.ts[["Q"]], type = "l", main = gauge, xlab = "Date", ylab = "Q")
            
            points(flow.ts[which(high.flows == 1), "Date"], flow.ts[which(high.flows == 1), "Q"], col = "red", cex = 0.25)
            
            abline(h = flow.threshold)
        }
        
    }
    
    if (ann.stats == F) {
        return(data.frame(high.spell.threshold = flow.threshold, n.events = n.events, spell.frequency = flood.frequency, ari = 1/flood.frequency, min.high.spell.duration = min.duration, avg.high.spell.duration = avg.duration, 
            med.high.spell.duration = med.duration, max.high.spell.duration = max.duration, avg.spell.volume = mean(spell.volumes, na.rm = T), avg.spell.peak = high.flow.av, sd.spell.peak = high.flow.sd, 
            avg.rise = mean.rise, avg.fall = mean.fall))
    } else {
        return(data.frame(high.spell.threshold = flow.threshold, n.events = n.events, spell.freq = flood.frequency, ari = 1/flood.frequency, min.high.spell.duration = min.duration, avg.high.spell.duration = avg.duration, 
            med.high.spell.duration = med.duration, max.high.spell.duration = max.duration, avg.spell.volume = mean(spell.volumes, na.rm = T), avg.spell.peak = high.flow.av, sd.spell.peak = high.flow.sd, 
            avg.rise = mean.rise, avg.fall = mean.fall, avg.max.ann = ann.maxs.mean, cv.max.ann = cv.max.ann, flood.skewness = flood.skewness, ann.max.timing = avg.max.day[["mean.doy"]], ann.max.timing.sd = avg.max.day[["sd.doy"]], 
            ann.max.min.dur = min.max.ann.duration, ann.max.avg.dur = avg.max.ann.duration, ann.max.max.dur = max.max.ann.duration, ann.max.cv.dur = cv.ann.duration))
    }
    
    
} 
