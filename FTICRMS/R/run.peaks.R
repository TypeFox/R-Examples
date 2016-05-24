`run.peaks` <-
function(trans.method = c("shiftedlog", "glog", "none"), add.par = 0, 
        subtract.base = FALSE, root.dir = ".", base.dir, peak.dir, 
        overwrite = FALSE, use.par.file = FALSE, par.file = "parameters.RData", 
        num.pts = 5, R2.thresh = 0.98, oneside.min = 1, 
        peak.method = c("parabola", "locmaxes"), calc.all.peaks = FALSE, 
        gengamma.quantiles = TRUE, peak.thresh = 3.798194){
    fail <- 0
    zeros <- 0
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baselines", sep="")}
    if(missing(peak.dir)){peak.dir <- paste(root.dir, "/All_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
        tmp <- match.call()
        tmp[[1]] <- as.name("list")
        tmp <- eval(tmp)
        if(length(tmp) > 0){
            for(i in 1:length(tmp)){
                assign(names(tmp)[i],tmp[[i]])
            }
        }
    }
    peak.method <- match.arg(peak.method)
    trans.method <- match.arg(trans.method)
    if(!file.exists(peak.dir)){
        dir.create(peak.dir)
    }
    for(j in list.files(base.dir)){
        if(!file.exists(paste(peak.dir, "/", sub("\\.RData$", "_peaks.RData", j), sep="")) || 
                overwrite){
            load(paste(base.dir, "/", j, sep=""))
            names(spect) <- c("Freq","LogAmp")
            if(subtract.base){
                spect$LogAmp <- spect$LogAmp-spect.base
            } else if(any(spect$LogAmp == 0)){
                zeros <- zeros + 1
                spect$LogAmp[spect$LogAmp == 0] <- min(spect$LogAmp[spect$LogAmp > 0])
            }
            if(trans.method=="shiftedlog"){
                if(subtract.base){
                    add.par <- add.par - min(spect$LogAmp)
                }
                spect$LogAmp <- log(spect$LogAmp+add.par)
            }
            if(trans.method=="glog"){
                spect$LogAmp <- log((spect$LogAmp+sqrt(add.par + spect$LogAmp^2))/2)
            }
            if(calc.all.peaks){
                thresh <- -Inf
                spect <- spect[order(spect$Freq),]
                numsplit <- (dim(spect)[1] %/% 100000) + 1
                splits <- (1:(numsplit-1)) * floor(dim(spect)[1]/numsplit)
                splits <- apply(cbind(c(1,splits-(num.pts-oneside.min-1)),
                    c(splits+(num.pts-oneside.min-1),dim(spect)[1])), 1, list)
                splits <- lapply(splits, function(x){seq(x[[1]][1],x[[1]][2])})
                for(i in 1:numsplit){
                    all.peaks.tmp <- locate.peaks(spect[splits[[i]],], num.pts=num.pts,
                        R2.thresh=R2.thresh, oneside.min=oneside.min,
                        peak.method=peak.method, thresh=thresh)
                    save(all.peaks.tmp, file=paste(peak.dir, "/", sub("\\.RData",
                        paste("_peaks", i, ".RData", sep=""), j), sep=""))
                    rm(all.peaks.tmp)
                }
                for(i in numsplit:1){
                    load(paste(peak.dir, "/", sub("\\.RData", paste("_peaks", i, ".RData",
                        sep=""), j), sep=""))
                    if(i == numsplit){
                        all.peaks <- all.peaks.tmp
                    } else {
                        all.peaks <- rbind(all.peaks.tmp,all.peaks)
                    }
                    rm(all.peaks.tmp)
                    file.remove(paste(peak.dir, "/", sub("\\.RData", paste("_peaks", i,
                        ".RData", sep=""), j), sep=""))
                }
            } else {
                if(gengamma.quantiles){
                    if(subtract.base){
                        spect.base <- 0
                    }
                    if(trans.method=="shiftedlog"){
                        thresh <- log(peak.thresh) + log(spect.base+add.par)
                    }
                    if(trans.method=="glog"){
                        thresh <- log(peak.thresh) +
                            log((spect.base+sqrt(add.par + spect.base^2))/2)
                    }
                } else {
                    thresh <- .biweight.FTICRMS(spect$LogAmp, K=3*peak.thresh/2)
                    thresh <- thresh$center + 3*peak.thresh/2*thresh$scale
                }
                all.peaks <- locate.peaks(spect, num.pts=num.pts,
                        R2.thresh=R2.thresh, oneside.min=oneside.min,
                        peak.method=peak.method, thresh=thresh)
            }
            all.peaks <- unique(all.peaks)
            rownames(all.peaks) <- 1:dim(all.peaks)[1]
            save(all.peaks, file=paste(peak.dir, "/", sub("\\.RData$",
                "_peaks.RData", j), sep=""))
        } else {
            fail <- fail + 1
        }
    }
    if(fail) {
        warning(paste(fail, "peak file(s) already existed and overwrite = FALSE; those file(s) not updated"))
    }
    if(zeros){
        warning(paste(zeros, ifelse(zeros==1, "spectrum", "spectra"), "had one or more zero entries for amplitude; those entries replaced by minimum value of amplitude"))
    }
}

