"importJDBC" <- 
function(sqlQuery,
         table,
         driverClass=sjdbcOptions()$driverClass,
         con=sjdbcOptions()$con, 
         user=sjdbcOptions()$user,
         password=sjdbcOptions()$password, 
         keepAlive=sjdbcOptions()$keepAlive,
         bigdata=FALSE)
{	
	if (bigdata == TRUE) stop("bigdata not supported in this version")
	if (driverClass == "") stop("driverClass is required.")
	if (con == "") stop("con is required.")
	if (missing(sqlQuery) && missing(table)) stop("Either sqlQuery or table is required.")
	if (missing(sqlQuery)) sqlQuery = paste("SELECT * FROM", table)

	# Specify the id for the Java hashtable.
	handleid <- "sqljdbcid"
	
	# Execute the query and point to the result in Java.
	.JavaMethod("SJDBCBridge", "sjdbcImportData", 
		"(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;)V",
		driverClass, con, user, password, sqlQuery, handleid)
		
	# make sure we unregister result set on exit
	on.exit(.JavaMethod("SJDBCResultSetUtilities", "unregister", "(Ljava/lang/String;)V",
					  handleid))


	# Retrieve the result as an data.frame.  
	res <- sjdbcGetResultSet(handleid, default.num.rows=100)

	if (!keepAlive) sjdbcCloseConnection()
	res
}