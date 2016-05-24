mdbExistsTable <- function(con, name){
	if(class(con)=="PostgreSQLConnection"){
		return(dbExistsTable(con, tolower(name)))
	} else {
		return(dbExistsTable(con, name))
	}
}
