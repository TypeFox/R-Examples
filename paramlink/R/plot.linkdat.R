
.as.kinship2_pedigree = function(x, deceased=numeric(0), aff2=NULL) {
    ped = relabel(x$pedigree, x$orig.ids) 
    AFF = ped[, 'AFF']; AFF[AFF == 0] = NA; AFF = AFF - 1 #kinship2:::plot.pedigree uses aff codes 1=aff, 0=unaff, NA=unknown
    if(!is.null(aff2)) {
        stopifnot(is.numeric(aff2), length(aff2)==x$nInd, all(aff2 %in% 0:2))
        aff2[aff2==0] = NA; aff2 = aff2 - 1
        AFF = cbind(aff1 = AFF, aff2 = aff2)
    }
    
    status = ifelse(x$orig.ids %in% deceased, 1, 0)
    pedigree(id=ped[,'ID'], dadid=ped[,'FID'], momid=ped[,'MID'], sex=ped[,'SEX'], affected=AFF, status=status)
}

plot.linkdat <- function(x, marker=NULL, alleles=NULL, sep="/", missing="-", skip.empty.genotypes=FALSE,
                         id.labels=x$orig.ids, title=paste("Family", x$famid), available=FALSE, col = 1,
                         deceased=numeric(0), starred=numeric(0), aff2=NULL, margins=c(0.6, 1, 4.1, 1), ...) {
    
    pedigree = .as.kinship2_pedigree(x, deceased=deceased, aff2=aff2)
    
    if(available) 
        cols = ifelse(x$orig.ids %in% x$available, 2, 1)
    else 
        cols = col
    
    if(!is.null(lb <- x$loop_breakers) && length(id.labels)==x$nInd) {origint = .internalID(x, lb[,1]); copyint = .internalID(x, lb[,2]); id.labels[copyint] = paste(id.labels[copyint], id.labels[origint], sep="=")}
    
    strid = rep(id.labels, length.out=x$nInd)
    starred = .internalID(x, starred)
    strid[starred] = paste(strid[starred], "*", sep="")
    if (!is.null(marker)) {
        if (inherits(marker, "marker")) 
            m = list(marker)
        else if(is.list(marker) && all(sapply(marker, inherits, what="marker"))) 
            m = marker  #marker must be either a 'marker' object, a list of such, or an integer vector.
        else 
            m = x$markerdata[marker]
        #chrom = ifelse(all(sapply(m, function(mm) identical(23L, as.integer(attr(mm, 'chrom'))))), 'X', 'AUTOSOMAL')
        gg = .prettyMarkers(m, alleles=alleles, sep=sep, missing=missing, singleCol=TRUE, sex=x$pedigree[, 'SEX'])
        geno = apply(gg, 1, paste, collapse="\n")
        if (skip.empty.genotypes) 
            geno[rowSums(do.call(cbind, m))==0] = ""
        if (is.null(strid) || !any(nzchar(strid))) 
            strid = geno 
        else 
            strid = paste(strid, geno, sep="\n")
    }    
    oldmar = par(mar=margins)  # without this, par() does not equal 'margins'...(why??) Needed for centered title.
    plot(pedigree, id=strid, col=cols, mar=margins, ...)
    if(!is.null(title)) 
        title(title)
    par(oldmar)
}


plot.singleton = function(x, marker=NULL, alleles=NULL, sep="/", missing="-", 
    id.labels=x$orig.ids, title=paste("Family", x$famid), available=FALSE, col = 1,
    starred=numeric(0), margins=c(4.1, 1, 4.1, 1), symbolsize=4, offset=0.5+0.2*symbolsize, textsize=1, ...) {
   
    sex=x$pedigree[, 'SEX']
   
    if(available) 
        cols = ifelse(x$orig.ids %in% x$available, 2, 1)
    else cols = col
   
    strid = rep(id.labels, length.out=x$nInd)
    starred = .internalID(x, starred)
    strid[starred] = paste(strid[starred], "*", sep="")
    if (!is.null(marker)) {
        if(inherits(marker, "marker")) 
            m = list(marker)
        else if(is.list(marker) && all(sapply(marker, inherits, what="marker")))
            m = marker  #marker must be either a 'marker' object, a list of such, or an integer vector.
        else 
            m = x$markerdata[marker]
        #chrom = ifelse(all(sapply(m, function(mm) identical(23L, as.integer(attr(mm, 'chrom'))))), 'X', 'AUTOSOMAL')
        gg = .prettyMarkers(m, alleles=alleles, sep=sep, missing=missing, singleCol=TRUE, sex=sex)
        geno = gg[,1]
        for(i in seq_along(m)[-1]) 
            geno = paste(geno, gg[, i], sep="\n")
        if(is.null(strid) || !any(nzchar(strid))) 
            strid = geno 
        else strid = paste(strid, geno, sep="\n")
    }
   
    plot(1, 1, pch=ifelse(sex==1, 0, 1), col=cols, axes=FALSE, xlab="", ylab="", cex=symbolsize, mar=margins, ...)
    text(1, 1, strid, pos=1, offset=offset, cex=textsize)
    if(!is.null(title)) title(title)
}

