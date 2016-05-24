get.stat.data <-
function(app.id, param, stat.id){
	stat.data.value.api.uri <- get.api.url(app.id, 'stat-data-value', param, stat.id)
	stat.data.note.api.uri  <- get.api.url(app.id, 'stat-data-note',  param, stat.id)
	value.list.names <- strsplit(readLines(stat.data.value.api.uri,1), "\t")[[1]]
	note.list.names  <- strsplit(readLines(stat.data.note.api.uri,1), "\t")[[1]]
	value.data.list <- rep("", length(value.list.names))
	note.data.list  <- rep("", length(note.list.names))
	names(value.data.list) <- value.list.names
	names(note.data.list)  <- note.list.names
	stat.data.value <- scan(stat.data.value.api.uri, as.list(value.data.list), skip=1, sep="\t")
	stat.data.value$value <- as.numeric(stat.data.value$value)
	stat.data.note <- scan(stat.data.note.api.uri, as.list(note.data.list), skip=1, sep="\t")
	stat.data <- list(
		value = as.data.frame(stat.data.value),
		note  = stat.data.note$value
	)
	return(stat.data)
}

