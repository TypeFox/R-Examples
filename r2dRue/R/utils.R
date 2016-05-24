# Authors: Robert J. Hijmans 
# International Rice Research Institute
# contact: r.hijmans@gmail.com
# Date : October 2008
# Version 0.9
# Licence GPL v3

readIniFile <- function(filename, token='=', commenttoken=';') {
	
	trim<- function(x){
		f <- function(s) {return( gsub('^[[:space:]]+', '',  gsub('[[:space:]]+$', '', s) ) )}
		return(unlist(lapply(x, f)))
	}

	strSplitOnFirstToken <- function(s, token="=") {
		pos <- which(strsplit(s, '')[[1]]==token)[1]
		if (is.na(pos)) {
			return(c(trim(s), NA)) 
		} else {
			first <- substr(s, 1, (pos-1))
			second <- substr(s, (pos+1), nchar(s))
			return(trim(c(first, second)))
		}
	}
	
	strsp <- function(s){strSplitOnFirstToken(s, token=token)}
	
	strSplitComment <- function(s,  token=";") { 
		# ";" is the start of a comment .
		strSplitOnFirstToken(s, token=";") 
	}
	strspcom <- function(s){strSplitComment(s, token=commenttoken)}
	
	
	if (!file.exists(filename)) { stop(paste(filename, " does not exist")) }
	
	Lines <- readLines(filename,  warn = FALSE)
	Lines <- trim(Lines)
	
	ini <- lapply(Lines, strspcom) 
	
	Lines <- matrix(unlist(ini), ncol=2, byrow=TRUE)[,1]
	ini <- lapply(Lines, strsp) 
	
	ini <- matrix(unlist(ini), ncol=2, byrow=TRUE)
	ini <- subset(ini, ini[,1] != "")
	
	ns <- length(which(is.na(ini[,2])))
	if (ns > 0) {
		sections <- c(which(is.na(ini[,2])), length(ini[,2]))
		
# here I should check whether the section text is enclused in [ ]. If not, it is junk text that should be removed, rather than used as a section
		ini <- cbind("", ini)
		for (i in 1:(length(sections)-1)) {
			ini[sections[i]:(sections[i+1]), 1] <- ini[sections[i],2]
		}	
		ini[,1] <- gsub("\\[", "", ini[,1])
		ini[,1] <- gsub("\\]", "", ini[,1])
		sections <- sections[1:(length(sections)-1)]
		ini <- ini[-sections,]
	} else {
		ini <- cbind("", ini)	
	}
	
	colnames(ini) <- c("section", "name", "value")
	########### como dataframe
	out=as.data.frame(t(ini[,3]),stringsAsFactors=FALSE)
	names(out)=ini[,2]
	return(out)
}

.gdalWriteFormats <- function() {
	gd <- gdalDrivers()
	gd <- as.matrix(subset(gd, gd[,3] == T))
	i <- which(gd[,1] %in% c('VRT', 'MEM', 'MFF', 'MFF2'))
	gd[-i,]
}


isSupportedGDALFormat <- function(dname) {
	if (!require(rgdal)) { stop() }
	gd <- .gdalWriteFormats()
	res <- dname %in% gd[,1]	
	return(res)
}

