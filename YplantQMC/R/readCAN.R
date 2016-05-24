# Reads .CAN format
readCAN <- function(canfile){

	r <- readLines(canfile)
	isf_row <- grep("INDIRECT SITE FACTOR",r)
	
	risf <- r[(isf_row+2):(isf_row+21)]
	# get rid of first 'column', some extra spaces.
	risf <- sapply(risf, function(x)substr(x, 6, nchar(x)))
	risf <- paste(risf, collapse="\n")
	hemimat <- read.table(text=risf)

	hemimat[,9] <- NULL
	
	# values are 'contributions to total'. Multiply by 160 to get gap fraction,
	# cf. 'plantcanopy' in Dassim.pas.
	hemimat <- hemimat * 160

	hemimat <- as.matrix(hemimat)
	colnames(hemimat) <- NULL
	
	# sort in ascending altitude.
	hemimat <- hemimat[nrow(hemimat):1,]
	
	# Relative size of sky sectors.
	weights <- diff(sin((pi/180)*seq(0,90,by=4.5)))
	weightm <- matrix(rep(weights, each=8), ncol=8, byrow=TRUE)
	
	# Openness relative to area of tile.
	hemimat <- hemimat/(20*weightm)
	
	# clip > 1 and < 0.
	hemimat[hemimat > 1] <- 1
	hemimat[hemimat < 0] <- 0
	
	
return(hemimat)
}