TraMineR.checkupdates <- function() {

  suppressWarnings(descr<-packageDescription("TraMineR"))
  if(!inherits(descr, "packageDescription")) stop(" [>] TraMineR is not installed")
  tmp <- tempfile(pattern = "file", tmpdir = tempdir())
  #suppressWarnings(z <- try(source("http://mephisto.unige.ch/traminer/VERSIONS", local=T), T))
  suppressWarnings(z <- try(download.file("http://mephisto.unige.ch/traminer/VERSIONS", destfile=tmp, quiet=T),T))
	
  if(inherits(z, "try-error")) {
	stop(" [>] Connection error")
  }
  
  else {
	arrayvers <- suppressWarnings(read.table(tmp, header=F))
	stable.chk <- as.character(arrayvers[1,1])
	dev.chk <- as.character(arrayvers[2,1])
	stable.vec <- extract.ver(stable.chk)
	dev.vec <- extract.ver(dev.chk)
	current.vec <- extract.ver(descr$Version)

	if(as.numeric(current.vec[2])%%2==0) {
		current.str <- "stable"
	}
	else current.str <- "development"

	## checks for the stable version
	if(current.vec[1] < stable.vec[1]) {
		message(paste(" [>] There is a new MAJOR stable version (", stable.chk, "), you have ", current.str, " version " , descr$Version, sep="")) 
	}
	else if(current.vec[2] < stable.vec[2]) {
		message(paste(" [>] There is a new MINOR stable version (", stable.chk, "), you have ", current.str, " version " , descr$Version, sep=""))
	}
	else if(length(stable.vec)==3) {
		if(length(current.vec)==3) {
			if(current.vec[3]<stable.vec[3]) {
				message(paste(" [>] There is a new bug fix for the stable version (", stable.chk, "), you have version " , descr$Version, sep=""))
			}
		}
		else if (current.vec[1]==stable.vec[1] & current.vec[2]==stable.vec[2]) {
			message(paste(" [>] There is a new bug fix for the stable version (", stable.chk, "), you have version " , descr$Version, sep=""))

		}
		else {
			cur.minor <- current.vec[2]
			if ((as.numeric(cur.minor)%%2)==0) message(paste(" [>] Your stable version", descr$Version, "is up-to-date."))
		}
	}
	


	## check for the dev version
	if(current.vec[1] < dev.vec[1]) {
		message(paste(" [>] There is a new MAJOR dev version (", dev.chk, "), you have ", current.str, " version " , descr$Version, sep="")) 
	}
	else if(current.vec[2] < dev.vec[2]) {
		message(paste(" [>] There is a new MINOR dev version (", dev.chk, "), you have ", current.str, " version " , descr$Version, sep=""))
	}
	else if(length(dev.vec)==3) {
		if(length(current.vec)==3) {
			if(current.vec[3]<dev.vec[3]) {
				message(paste(" [>] There is a new bug fix for the dev version (", dev.chk, "), you have version " , descr$Version, sep=""))
			}
		}
		else if (current.vec[1]==dev.vec[1] & current.vec[2]==dev.vec[2]) {
			message(paste(" [>] There is a new bug fix for the dev version (", dev.chk, "), you have version " , descr$Version, sep=""))

		}
		else message(paste(" [>] Your development version", descr$Version, "is up-to-date."))
	}

	else {
		cur.minor <- current.vec[2]
		if ((as.numeric(cur.minor)%%2)!=0)	message(paste(" [>] Your development version", descr$Version, "is up-to-date."))
	}
  } 
}

extract.ver <- function(x) {

	ver.list <- strsplit(x, "\\.")
	ver.unit <- ver.list[[1]][1]
	if(length(grep("-", ver.list[[1]][2]))>0) {
		ver.dec.list <- strsplit(ver.list[[1]][2], "-")
		ver.dec <- ver.dec.list[[1]][1]
		ver.bug <- ver.dec.list[[1]][2]
		return(c(ver.unit, ver.dec, ver.bug))
	}

	ver.dec <- ver.list[[1]][2]
	return(c(ver.unit, ver.dec))
}