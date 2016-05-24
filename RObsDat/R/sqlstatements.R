sqlstatements <- function(object, term){
	if(term=="now"){
		backend <- class(object@con)
		if(backend=="SQLiteConnection"){
			datestring =  'datetime("now")'
		} else if(backend=="MySQLConnection") {
			datestring = "NOW()"
		} else if(backend=="PostgreSQLConnection") {
			datestring = "now()"
		} else {
			todo(paste("Implement backend ", backend))
			browser()
		}
	} else if(term=="last_id"){
		backend <- class(object@con)
		if(backend=="SQLiteConnection"){
			datestring =  'last_insert_rowid()'
		} else if(backend=="MySQLConnection") {
			datestring = "LAST_INSERT_ID()"
		} else if(backend=="PostgreSQLConnection") {
			todo(paste("Implement lastval backend ", backend))
			browser()
			datestring = "lastval()"
		} else {
			todo(paste("Implement backend ", backend))
			browser()
		}
	} else {
			todo(paste("Implement term ", term))
			browser()
	}
	return(datestring)
}

