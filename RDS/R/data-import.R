





#' Import data from the 'RDSAT' format as an \code{rds.data.frame}
#' @description This function imports RDSAT data files as \code{rds.data.frame} objects.
#' @param file the name of the file which the data are to be read from.
#'           If it
#'           does not contain an _absolute_ path, the file name is
#'           _relative_ to the current working directory, 'getwd()'.
#'           Tilde-expansion is performed where supported.  As from R
#'           2.10.0 this can be a compressed file (see 'file')
#' @param delim The seperator defining columns. <auto> will guess the
#' 		delimitor based on the file.
#' @param N The population size (Optional).
#' @examples 
#' fn <- paste0(path.package("RDS"),"/extdata/nyjazz.rdsat")
#' rd <- read.rdsat(fn)
#' plot(rd)
#' @export
read.rdsat <- function(file,delim=c("<auto>","\t"," ",","),N=NULL){
	
	# The only part that we need to parse are the first two lines, subsequent lines
	# should be a valid tab-delimited table.
	header <- readLines(file,n=2)	
	line.one <- header[1]
	
	if(length(grep(pattern="RDS",line.one)) == 0 & (length(grep(pattern="rds",line.one)) == 0)){
		warning("RDSAT file does not begin with \"RDS\" or \"rds\".")
	}	
	
	delim <- delim[1]
	if(delim == "<auto>"){
		line.two1 <- strsplit(header[2],split="\t")[[1]] 
		line.two2 <- strsplit(header[2],split=" ")[[1]] 
		line.two3 <- strsplit(header[2],split=",")[[1]]
		if(length(line.two1)>length(line.two2) && length(line.two1)> length(line.two3))
			delim <- "\t"
		else if(length(line.two2)>length(line.two3))
			delim <- " "
		else
			delim <- ","
	}
	line.two <-strsplit(header[2],split=delim)[[1]]	
	header.number.of.rows <- line.two[1]
	header.max.coupons <- as.integer(line.two[2])
	header.na.symbol <- line.two[3]
	header.column.names <- line.two[4:length(line.two)]
	header.named.columns <- header.column.names[header.column.names != ""]
	
	# The number of unnamed columns should be determined as follows:
	#
	#      First column is an id
	#      Second column is network size
	#      Third column is coupon received.
	#      Fourth through [Four + (header.max.coupons -1)]^th are coupons given.
	
	initial.column.names <- c("id","network.size","own.coupon",
			paste("coupon",1:header.max.coupons,sep="."))
	
	coupon.columns <- 4:(4 + header.max.coupons - 1) 
	
	if(delim==" ") delim=""
	rds.data <- utils::read.table(file,skip=2,sep=delim,na.strings=header.na.symbol,as.is=FALSE,stringsAsFactors=FALSE,nrows=1)
	ccc <- rep(NA,ncol(rds.data))
	ccc[c(1,3,coupon.columns)] <- "character"
	ccc[2] <- "numeric"
	rds.data <- utils::read.table(file,skip=2,sep=delim,na.strings=header.na.symbol,as.is=FALSE,stringsAsFactors=FALSE,colClasses=ccc)
	
	# Strip out most blank lines
	rds.data <- rds.data[!apply(is.na(rds.data[,1:3]),1,all),]
	
	colnames(rds.data) <- make.names(c(initial.column.names, header.named.columns),
			unique=TRUE)
	
	#convert to rds.data.frame
	rds.data$recruiter.id <- rid.from.coupons(rds.data, subject.coupon='own.coupon', 
			coupon.variables=paste("coupon",1:header.max.coupons,sep="."),subject.id="id")
	rds.data <- as.rds.data.frame(rds.data, id='id', recruiter.id='recruiter.id', 
			network.size='network.size', max.coupons=header.max.coupons,
			population.size=N)
	rds.data$seed <- get.seed.id(rds.data)
	rds.data$wave <- get.wave(rds.data)
	
	rds.data
}

#' Import data saved using write.rdsobj
#' @param file the name of the file which the data are to be read from.
#'           If it
#'           does not contain an _absolute_ path, the file name is
#'           _relative_ to the current working directory, 'getwd()'.
#'           Tilde-expansion is performed where supported.  As from R
#'           2.10.0 this can be a compressed file (see 'file')
#' @export
read.rdsobj <- function(file){
	obj <- try(readRDS(file), silent=TRUE)
	if(inherits(obj,"try-error")){obj <- dget(file)}
	if(!inherits(obj,"rds.data.frame"))
		stop("Object stored in file is not an rds.data.frame object")
	assert.valid.rds.data.frame(obj)
	obj
}



