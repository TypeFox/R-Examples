#rdb
byeShort <- function(X, percentReqd=0.80, Expected=288, ToCount="doy", TruncToCount=TRUE, By=c("year","doy")){
	#there is almost certainly a better way to do this
	#in addition to checking for missing values, also checks for duplicates (checks duplicates first, removes those, then checks for missing)
	# NOTE: see notes in addNAs about the presence of "doy" and "year" columns in the input to byeShort
	dups <- function(x) x[!duplicated(round(x[,ToCount],9)),]
	X <- ddply(X, setdiff(By, ToCount), dups)
	ByInd <- data.frame(X[,By], "IND"=1:nrow(X))
	which_nrow <- function(x){ c("Size"=nrow(x), "Start"=min(x[,"IND"]), "Stop"=max(x[,"IND"]))}
	Sizes <- ddply(trunc(ByInd), By, which_nrow)
	TooShort <- Sizes[which(Sizes[,"Size"] < Expected*percentReqd), c("Start","Stop")]
	Start2Stop <- function(x) x[1]:x[2] #these last two steps could probably be combined into an is.element() approach that would be simpler
	WaveTo <- unlist(apply(TooShort, MARGIN=1, FUN=Start2Stop), use.names=FALSE)
	print(paste("Points removed due to incomplete day or duplicated time step:",length(WaveTo), sep=" "))
	flush.console()
	if(length(WaveTo)!=0){
		Xr <- X[-WaveTo,]
	}else{
		Xr <- X
	}
	return(Xr)
}
