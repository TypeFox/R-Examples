"executeJDBC" <- 
function(sqlQuery, driverClass=sjdbcOptions()$driverClass, con=sjdbcOptions()$con, 
	 user=sjdbcOptions()$user, password=sjdbcOptions()$password, 
	 keepAlive=sjdbcOptions()$keepAlive)
{	
	if (driverClass == "") stop("driverClass is required.")
	if (con == "") stop("con is required.")
	if (missing(sqlQuery)) stop("sqlQuery is required.")
	
	# Execute the query
	res <- .JavaMethod("SJDBCBridge", "sjdbcExecute", 
		"(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;)I",
		driverClass, con, user, password, sqlQuery)
		
	if (!keepAlive) sjdbcCloseConnection()
	
	res
}
