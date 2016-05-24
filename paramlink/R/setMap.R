setMap = function(x, map, dat, pos=NULL, verbose=TRUE) {
	if(x$nMark==0) {x$map=NULL; return(x)}
	if(!is.null(pos)) {
		if(length(pos) != x$nMark) stop("Length of 'pos' argument does not match the number of markers.")
		x$map = data.frame(CHR=1, MARKER=paste("M", 1:x$nMark, sep=""), POS=pos, stringsAsFactors=FALSE)
		return(x)
	}
	
	stopifnot(!missing(map), !is.null(map), (is.data.frame(map) || is.character(map)))
	if(is.data.frame(map)) {
		if(ncol(map)>=3 && nrow(map)==x$nMark) {names(map)[1:3] = c('CHR', 'MARKER', 'POS'); x$map = map}
		else warning("Map not set: Something is wrong with the 'map' data frame.")
		return(x)
	}
	
	stopifnot(!missing(dat), !is.null(dat), is.character(dat))
	rawmap = read.table(map, as.is=TRUE, header=FALSE)
	if(!any(1:9 %in% strsplit(rawmap[1,1],"")[[1]])) rawmap = rawmap[-1,] #If no number occurs in first entry, first row is assumed to be header. 
	rawmap[[1]][rawmap[[1]]=="X"] = 23
	map1 = data.frame(CHR = as.numeric(rawmap[,1]), MARKER = as.character(rawmap[, 2]), POS = as.numeric(rawmap[, 3]), stringsAsFactors=FALSE) 
	rawdat = read.table(dat, as.is=TRUE, header=FALSE)
	dat = as.character(rawdat[rawdat[1]=="M", 2])  #names of all markers in map
		
	Mmatch = match(dat, map1$MARKER, nomatch=0)
	if(any(Mmatch == 0)) {
		del = dat[Mmatch==0]
		if(verbose) cat("Deleting the following marker(s), which are not found in the map file:\n", paste(del, collapse="\n"), "\n")
		x$markerdata[Mmatch==0] = NULL
	}

	map = map1[Mmatch, ]
	map = map[order(map$CHR, map$POS), ] 
		
	ord = match(dat, map$MARKER, nomatch=0)
	x$markerdata = x$markerdata[ord]
	x$nMark = length(x$markerdata)
	
	x$map = map
	x
}

readDatfile = function(datfile, chrom, write_to=NULL) {
	dat = lapply(strsplit(readLines(datfile), split=" |\t"), function(v) v[v!=""])
    #NB: comment.char = "<" necessary?
    nMark = as.numeric(dat[[1]])[1] - 1
	xlinked = as.numeric(dat[[1]])[3]
    ordering = as.numeric(dat[[3]])
    markernames = sapply(dat[6 + xlinked + (1:nMark)*2], '[', 4)
	pos = cumsum(as.numeric(dat[[length(dat)-1]]))
    if(ordering[1]==2) ordering = c(1,ordering)
    if(length(markernames) - length(pos) == 1) pos = c(0, pos)
    stopifnot(all(ordering==seq_len(nMark+1)), length(dat)==7 + xlinked + 2*nMark + 3)
	
	equal = (pos[-1]==pos[-length(pos)]); k=0
	for(i in 2:length(pos)) 
		if (equal[i-1]) pos[i]=pos[i] + 0.0001*(k <- k+1) else k=0  #if consecutive entries are equal, add 0.0001's. 
		
	map = data.frame(CHR=chrom, MARKER=markernames, POS=pos, stringsAsFactors=F)
	
	freqlist = lapply(dat[7 + xlinked + (1:nMark)*2], function(r) as.numeric(r))
    nalls = lengths(freqlist, use.names=F)
    L = sum(nalls) + length(nalls)
    cum = cumsum(c(1, nalls+1))
    length(cum) = length(nalls) #remove last
    col1 = rep("A", L)
    col1[cum] = "M"

    col2 = character(L)
    col2[cum] = markernames
    allalleles = unlist(lapply(nalls, seq_len))
    col2[-cum] = allalleles
    
    col3 = character(L)
    allfreqs = unlist(freqlist)
    col3[-cum] = format(allfreqs, scientifit=F, digits=6)
      
    freq = cbind(col1, col2, col3, deparse.level=0)
    
    merlindat = cbind(c("A", rep("M", nMark)), c("my_disease", markernames))
	if(!is.null(write_to)) {
        write.table(map, file=paste(write_to, "map", sep="."), row.names=F, col.names=F, quote=F)
        write.table(merlindat, file=paste(write_to, "dat", sep="."), row.names=F, col.names=F, quote=F)
        write.table(freq, file=paste(write_to, "freq", sep="."), row.names=F, col.names=F, quote=F)
    }   
	invisible(list(dat=merlindat, map=map, freq=freq))
}

.SNPfreq = function(markernames, Bfreq, allele1="1", allele2="2", file=NULL) {
    # Create and write freq file in extended Merlin format.
    # Bfreq numerical vector with same length as markernames (contains freqs for allele 2)
    
    stopifnot(length(markernames)==length(Bfreq))
    n = length(markernames)
    col1 = rep(c('M','A','A'), n)
    col2 = as.character(rbind(markernames, allele1, allele2))
    col3 = as.character(rbind("", 1-Bfreq, Bfreq))
    res = cbind(col1, col2, col3)
    if (!is.null(file)) write.table(res, file=file, col.names=F, row.names=F, quote=F)
    invisible(res)
}
