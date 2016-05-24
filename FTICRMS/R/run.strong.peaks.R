`run.strong.peaks` <-
function(cor.thresh = 0.8, isotope.dist = 7, pre.align = FALSE, 
        align.method = c("PL", "spline", "affine", "none"), align.fcn = NA,
        root.dir = ".", lrg.dir, lrg.file = "lrg_peaks.RData",
        overwrite = FALSE, use.par.file = FALSE, par.file = "parameters.RData"){
    if(missing(lrg.dir)){lrg.dir <- paste(root.dir, "/Large_Peaks", sep="")}
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
    align.method <- match.arg(align.method)
    load(paste(lrg.dir, "/", lrg.file, sep=""))
    if(!exists("amps") || overwrite){
        if(!identical(pre.align, FALSE)){
            if(!identical(class(pre.align),"list")){
                pre.align <- list(targets=median(lrg.peaks$Center_hat),
                    actual=data.frame(median(lrg.peaks$Center_hat)-pre.align),
                    align.method=align.method)
                if(use.par.file){
                    tmp <- extract.pars(par.file, root.dir)
                    tmp$pre.align <- pre.align
                    do.call(make.par.file, tmp)
                    load(paste(root.dir, "/", par.file, sep=""))
                }
            }
            if(!is.na(align.fcn)){
                pre.align$targets <- align.fcn$fcn(pre.align$targets)
                pre.align$actual <- align.fcn$fcn(pre.align$actual)
                lrg.peaks$Center_hat <- align.fcn$fcn(lrg.peaks$Center_hat)
            }
            pre.spl <- do.call(.alignment, pre.align)
            for(i in levels(lrg.peaks$File)){
                if(identical(class(pre.spl[[i]]),"lm")){
                    lrg.peaks$Center_hat[lrg.peaks$File==i] <- predict(pre.spl[[i]], 
                        data.frame(X=lrg.peaks$Center_hat[lrg.peaks$File==i]))
                } else {
                    lrg.peaks$Center_hat[lrg.peaks$File==i] <- predict(pre.spl[[i]], 
                        lrg.peaks$Center_hat[lrg.peaks$File==i])$y
                }
            }
            if(!is.na(align.fcn)){
                lrg.peaks$Center_hat <- align.fcn$inv(lrg.peaks$Center_hat)
            }
        }
        lrg.peaks <- lrg.peaks[order(lrg.peaks$Center_hat),]
        numfiles <- length(levels(lrg.peaks$File))
        strong <- which(sapply(1:(dim(lrg.peaks)[1]-numfiles), function(x){
            length(unique(lrg.peaks$File[x:(x+numfiles-1)])) == numfiles}))
        if(length(strong)==0){
            warning("No peaks were in all samples")
            amps <- NA
            centers <- NA
        } else {
            centers <- lapply(strong, function(x){
                lrg.peaks[x:(x+numfiles-1),c("Center_hat","File")]})
            centers <- data.frame(lapply(centers, function(x){x$Center_hat[order(x$File)]}))
            colnames(centers) <- sapply(centers,function(x){round(mean(x),3)})
            rownames(centers) <- levels(lrg.peaks$File)
            amps <- lapply(strong, function(x){
                lrg.peaks[x:(x+numfiles-1),c("Max_hat","File")]})
            amps <- data.frame(lapply(amps, function(x){x$Max_hat[order(x$File)]}))
            colnames(amps) <- colnames(centers)
            rownames(amps) <- levels(lrg.peaks$File)
            if(length(strong)>1){
                keep <- c(1,which(diff(strong)!=1)+1, dim(centers)[2]+1)
                for(i in 1:(length(keep)-1)){
                    keep[i] <- keep[i]-1+which.min(apply(centers[,keep[i]:(keep[i+1]-1),drop=FALSE],
                        2,function(x){diff(range(x))}))
                }
                keep <- keep[-length(keep)]
                centers <- centers[,keep,drop=FALSE]
                amps <- amps[,keep,drop=FALSE]
                amps <- amps[,sapply(centers,function(x)diff(range(x)))<.1,drop=FALSE]
                centers <- centers[,sapply(centers,function(x)diff(range(x)))<.1,drop=FALSE]
                if(dim(amps)[2]==0){
                    warning("No peaks were in all samples")
                    amps <- NA
                    centers <- NA
                } else {
                    cors <- cor(amps)
                    isotopes <-  outer(sapply(centers,mean),sapply(centers,mean), "-")
                    isotopes <- cors>=cor.thresh & abs(isotopes)<isotope.dist
                    isotopes <- apply(isotopes & upper.tri(isotopes),2,any)
                    centers <- centers[,!isotopes,drop=FALSE]
                    amps <- amps[,!isotopes,drop=FALSE]
                }
            }
        }
        save(lrg.peaks, amps, centers, file=paste(lrg.dir, "/", lrg.file, sep=""))
    } else {
        warning("File created by previous run of run.strong.peaks exists and overwrite = FALSE; file not updated")
    }
}
