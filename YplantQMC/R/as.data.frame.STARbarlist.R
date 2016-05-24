#'@method as.data.frame STARbarlist
#'@S3method as.data.frame STARbarlist
as.data.frame.STARbarlist <- function(x, ...){
	res <- lapply(x, function(x)unlist(x[1:5]))
	m <- do.call("rbind",res)
	dfr <- as.data.frame(m)
	dfr$pfile <- sapply(x, "[[", "pfile")
	dfr$lfile <- sapply(x, "[[", "lfile")

	return(dfr)
}