#'Reads the data from a .p file.
#'
#'Reads all data of the plant input file (\code{.p}) into a dataframe.
#'
#'
#'@param pfile Name of the p file (a character string).
#'@return Dataframe. The column names correspond to the labels in the .p file,
#'but duplicates have a numeric subscript;
#'
#'\describe{ \item{N}{ - node number} \item{MN}{ - mother node} \item{ste}{ -
#'stem (1) or branch (2)} \item{Lt}{ - leaf type} \item{Az, Az.1, Az.2}{ -
#'Azimuth of stem, branch, petiole and leaf.} \item{An, An.1, An.2}{ - Angle of
#'stem, branch, petiole and leaf.} \item{L, L.1, L.2, L.3}{ - Length of stem,
#'branch, petiole and leaf.} } More detail on the leaf angles:
#'
#'\describe{ \item{Az.3}{ - Leaf azimuth (azimuth angle of the normal to the
#'leaf surface).} \item{Or}{ - Leaf orientation (azimuth angle of the midrib).}
#'\item{An.3}{ - Leaf angle (angle between normal to the leaf surface and the
#'z-axis).}
#'
#'}
#'@author Remko Duursma
#'@keywords misc
#'@export
readp <- function(pfile=NA){

	if(is.na(pfile)){
		if(.Platform$OS.type != "windows" || !interactive())
			stop("Please provide a plant (.P) file")
		pfile <- file.choose()
	}

    r <- tolower(readLines(pfile, warn=FALSE))
    r <- gsub("[[:space:]]", "", r)
	
	# Find labels (N MN STE ...)
	headerline <- grep("nmnste",r)
	
	# If still no success, find the line that starts with 0 and 0
	if(any(substr(r, 1,2) == "00")){
		headerline <- which(substr(r, 1,2) == "00")[1] - 1
	}
	
	# If not found, find previous line that says "Petiole", amongst other things
	if(length(headerline) == 0){
		gr <- grep("petiole", r)
		if(length(gr) > 0)headerline <- gr[length(gr)] + 1
	}
	
	# If still no success, assume header is in 5th line (kind of defaultish location)
	if(length(headerline)==0)headerline <- 5
	
	pdata <- read.table(pfile, skip=headerline)
	
	# If not read correctly, try this;
	if(ncol(pdata) != 20)
		 pdata <- read.table(pfile, skip=headerline-1, fill=TRUE,header=TRUE)
	
	if(ncol(pdata) != 20)stop("Need 20 columns exactly.\n")
	
	pdata[is.na(pdata)] <- 0
	
	names(pdata) <- c("N", "MN", "ste", "Lt", "Az", "An", "L", "D", "Az.1", "An.1", "L.1",
                      "D.1", "Az.2", "An.2", "L.2", "D.2", "Or", "Az.3", "An.3", "L.3")

pdata
}
