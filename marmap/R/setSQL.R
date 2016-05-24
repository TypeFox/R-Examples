setSQL = function(bathy, header = TRUE, sep = ",", db.name = "bathy_db"){

	# prepare ("connect") SQL database
	con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db.name)
	# data frame -> database table.
	DBI::dbWriteTable(con, name="bathy_db", value=bathy, header=header, sep=sep)

}