subsetSQL = function(min_lon, max_lon, min_lat, max_lat, db.name="bathy_db"){

	# prepare ("connect") SQL database
	con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db.name)

	cn <- DBI::dbListFields(con,"bathy_db")

	# build SQL request: 
    paste("SELECT * from bathy_db where ", cn[1], " >", min_lon,
    	  "and ", cn[1], " <", max_lon," and ", cn[2], " >",min_lat," and ", cn[2], " <",max_lat) -> REQUEST
    
    # send request and retrieve results
    res <- DBI::dbSendQuery(con, REQUEST)
    data <- DBI::fetch(res, n = -1)
	DBI::dbClearResult(res)
	DBI::dbDisconnect(con)
	if (!is.numeric(data[,3])) data[,3] <- suppressWarnings(as.numeric(data[,3]))
    return(as.bathy(data))
}