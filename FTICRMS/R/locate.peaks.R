`locate.peaks` <-
function(peak.base, num.pts = 5, R2.thresh = 0.98, oneside.min = 1,
        peak.method = c("parabola","locmaxes"), thresh = -Inf){
    loc.max <- .loc.maxes(peak.base[,2])
    loc.max <- intersect(which(peak.base[,2]>=thresh), loc.max)
    peak.method <- match.arg(peak.method)
    if(peak.method=="parabola"){
        loc.max <- loc.max[loc.max >= oneside.min+1 & loc.max <= length(peak.base[,2])-oneside.min]
        all.peaks <- lapply(loc.max,function(x){
            locs <- x + ((-num.pts+oneside.min+1):(num.pts-oneside.min-1))
            locs <- locs[locs>=1 & locs<=dim(peak.base)[1]]
            .peak.parab(peak.base[locs,], num.pts, R2.thresh)
        })
        if(length(all.peaks)>0){
            all.peaks <- data.frame(do.call(rbind,all.peaks))
            all.peaks <- all.peaks[all.peaks[,3]>0,,drop=F]
        } else {
            all.peaks <- data.frame(matrix(0,nrow=0,ncol=3))    
        }
        names(all.peaks) <- c("Center_hat","Max_hat","Width_hat")     
        if(dim(all.peaks)[1]>0){
            rownames(all.peaks) <- 1:dim(all.peaks)[1]
            for(i in 1:3){
                all.peaks[,i] <- as.numeric(all.peaks[,i])
            }
        }
    } else if(peak.method == "locmaxes"){
        if(length(loc.max)>0){
            all.peaks <- data.frame(peak.base[loc.max,], NA)
            rownames(all.peaks) <- 1:dim(all.peaks)[1]
        } else {
            all.peaks <- data.frame(matrix(0,nrow=0,ncol=3))
        }
        names(all.peaks) <- c("Center_hat","Max_hat","Width_hat")     
    }
    all.peaks
}
