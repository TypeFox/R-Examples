
# Obsolete! Reads the E() module; but now all interesting stuff is in the qmcarea.out file.
readQMCout <- function(outfile){
	options(warn= -1)
	on.exit(options(warn=0))
  
	r <- readLines(outfile)
	n <- length(r)
	
	# g <- gsub("[^[:digit:].]+"," ",r)
	# linereads <- sapply(g, function(x)scan(textConnection(x), quiet=TRUE))
	# nnum <- sapply(linereads, length)
	
	# Remove everything that is not within the parentheses...
	r2 <- gsub(".+\\(", "", r)
	r2 <- gsub("\\).+", "", r2)
	r2 <- gsub("\\)", "", r2)
	
	pfd <- as.numeric(r2)

	id <- gsub(" E\\(.+","",r)
	id <- gsub("\001","",id)
	id <- as.numeric(id)
	
# return only vector of PFDs, for now. If things go bad, return ID as well to check...
return(pfd[!is.na(pfd)])
}
