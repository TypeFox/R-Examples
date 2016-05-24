char2seed <- function(x,set=TRUE,...){

	tmp <- c(0:9,0:25,0:25)
	names(tmp) <- c(0:9,letters,LETTERS)

	x <- gsub("[^0-9a-zA-Z]","",as.character(x))

	xsplit <- tmp[ strsplit(x,'')[[1]] ]

	seed <- sum(rev( 7^(seq(along=xsplit)-1) ) * xsplit)
        seed <- as.integer( seed %% (2^31-1) )

	if(set){
		set.seed(seed,...)
		return(invisible(seed))
	} else {
		return(seed)
	}
}

