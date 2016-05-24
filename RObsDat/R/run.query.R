run.query <- function(object, query){
	#replace " by ' because this is required by postgres
	queryA <- gsub("\"", "'", query)
	if(getOption("verbose.queries", default=FALSE)) print(queryA)

	IdbState(object)
	res <- dbGetQuery(object@con, queryA)
	return(res)
}
