# summary for maker maps

summaryGenMap <- function(map){

     # information from arguments
     if(class(map)=="gpData") map <- map$map
     if(is.null(map)) stop("No map available")
     chr <- map$chr
     pos <- map$pos

     # extract information
     # number of markers
     len <- tapply(pos,chr,length)
     rge <- tapply(pos,chr,max,na.rm=TRUE)-tapply(pos,chr,min,na.rm=TRUE)

     # differences of markers on each chr
     diffs <- tapply(pos[!is.na(pos)],chr[!is.na(pos)],diff,na.rm=TRUE)
     avDist <- as.numeric(lapply(diffs,mean,na.rm=TRUE))
     maxDist <- as.numeric(lapply(diffs,max,na.rm=TRUE))
     minDist <- as.numeric(lapply(diffs,min,na.rm=TRUE))

     # return data.frame
     ret <- data.frame(noM = len, length=rge, avDist = avDist, maxDist = maxDist, minDist = minDist,row.names=names(len))
     # keep same order as original map
     ret <- ret[order(order(unique(chr))),]
     
     # sum over all chr
     all <- data.frame(noM = sum(ret$noM,na.rm=TRUE),length= sum(ret$length,na.rm=TRUE),avDist=weighted.mean(ret$avDist,ret$noM,na.rm=TRUE),maxDist=max(ret$maxDist,na.rm=TRUE),minDist=min(ret$minDist,na.rm=TRUE))
     rownames(all) <- paste(rownames(ret)[1],"-",rownames(ret)[nrow(ret)])
     
     ret <- rbind(ret,all)
     
     return(ret)

}
