get.stat.list <-
function(app.id, param){
	stat.list.api.uri <- get.api.url(app.id, stat.type = 'stat-list', param)
	stat.list.names <- strsplit(readLines(stat.list.api.uri,1), "\t")[[1]]
	stat.list.value  <- rep("", length(stat.list.names))
	names(stat.list.value) <- stat.list.names
	stat.list <- scan(stat.list.api.uri, as.list(stat.list.value), skip=1, sep="\t")
	return(as.data.frame(stat.list))
}

