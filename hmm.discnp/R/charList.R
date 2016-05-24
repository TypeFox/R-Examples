charList <- function(y) {
if(is.matrix(y)) {
	warning(paste("Presenting \"y\" as a matrix is deprecated.\n",
                      "Change to presenting \"y\" either as a vector\n",
                      "or a list of vectors.\n"))
	y <- as.list(as.data.frame(y))
# Remove NA padding at the ends of the vectors in y.
	y <- lapply(y,function(x){
			if(all(is.na(x))) return(NULL)
			m <- max(which(!is.na(x)))
			x[1:m]
		      })
}
if(!is.list(y)) y <- list(y)
uval <- as.character(sort(unique(unlist(y))))
y    <- lapply(y,as.character)
if("missing" %in% uval) {
    whinge <- paste("At least one observation consists of the character\n",
                    "string \"missing\". This is a *reserved* word in the\n",
                    "hmm.discnp package and may not be used as an\n",
                    "observation value.\n")
    stop(whinge)
}

attr(y,"uval") <- uval
y
}
