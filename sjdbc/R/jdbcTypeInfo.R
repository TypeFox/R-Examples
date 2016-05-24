"jdbcTypeInfo"<-
function(driverClass = sjdbcOptions()$driverClass, con = 
	sjdbcOptions()$con, user = sjdbcOptions()$user, password = sjdbcOptions(
	)$password, keepAlive = sjdbcOptions()$keepAlive)
{
	if(driverClass == "")
		stop("driverClass is required.")
	if(con == "")
		stop("con is required.")
	# Specify the id for the Java hashtable.
	handleid <- "sqljdbcid"
	# Execute the query and point to the result in Java.
	.JavaMethod("SJDBCBridge", "sjdbcTypeInfo", 
		"(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;)V",
		driverClass, con, user, password, handleid)
	# Retreive the result as an data.frame.  
	res <- sjdbcGetResultSet(handleid)
	if(!keepAlive)
		sjdbcCloseConnection()
	res
}
