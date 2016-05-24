## Returns the longest suffix of path in a PST

lsuffix <- function(PST, path, asc=FALSE) {

	if (length(path)==1) { path <- seqdecomp(path) }
	sl <- length(path)
	
	if (asc) {
		lsuf <- sl
		while (paste(path[lsuf:sl], collapse="-") %in% rownames(PST[[sl-lsuf+1]]) & lsuf > 1) {
			lsuf <- lsuf-1
		}
	} else {
		lsuf <- 1
		while (!paste(path[lsuf:sl], collapse="-") %in% rownames(PST[[sl-lsuf+1]] ) & lsuf < sl) {
			lsuf <- lsuf+1
		}
	}

	message(" [>] longest suffix in the tree is: ", paste(rev(path[1:lsuf]), collapse="-"))

	return(path[1:lsuf])
}


