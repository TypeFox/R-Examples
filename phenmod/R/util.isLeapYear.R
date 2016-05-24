util.isLeapYear <- function(year){
	res <- vector(mode="logical", length=length(year))
	res <- rep(FALSE, length(res))
	res[which((year %% 4) == 0)] <- TRUE
	res[which((year %% 100) == 0)] <- FALSE
	res[which((year %% 400) == 0)] <- TRUE
	return(res)
}