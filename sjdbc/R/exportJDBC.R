"exportJDBC" <- 
function(data, table, appendToTable=TRUE, 
	 driverClass=sjdbcOptions()$driverClass, con=sjdbcOptions()$con, 
	 user=sjdbcOptions()$user, password=sjdbcOptions()$password, 
	 keepAlive=sjdbcOptions()$keepAlive, preserveColumnCase=FALSE,
	 batchSize = sjdbcOptions()$batchSize) {
	 
	if (driverClass == "") stop("driverClass is required.")
	if (con == "") stop("con is required.")
	if (missing(data)) stop("data is required.")
	if (missing(table)) stop("table is required.")	

	dfid <- table
	
  send.data.to.sjdbc <- function(dfid, data) {
    if (is.null(data)) {
      .JavaMethod("SJDBCBridge", "sjdbcClearData", "(Ljava/lang/String;)V", dfid)
      return()
    }
    .JavaMethod("SJDBCBridge", "sjdbcCreateData", "(Ljava/lang/String;)V", dfid)
    for (i in 1:length(data)) {
      if (inherits(data[,i], "integer")) {
		.JavaMethod("SJDBCBridge", "sjdbcAddIntegerColumn", "(Ljava/lang/String;Ljava/lang/String;[I)V",
				dfid, names(data)[i], data[,i])
      } else if (inherits(data[,i], "numeric")) {
        # numeric type
        .JavaMethod("SJDBCBridge", "sjdbcAddDoubleColumn", "(Ljava/lang/String;Ljava/lang/String;[D)V",
                    dfid, names(data)[i], data[,i])
      } else if (inherits(data[,i], c("Date", "POSIXlt", "POSIXct"))) {
        # numeric type
        .JavaMethod("SJDBCBridge", "sjdbcAddDateColumn", "(Ljava/lang/String;Ljava/lang/String;[Ljava/lang/String;)V",
                    dfid, names(data)[i], jdbcTimeDate(data[,i]))	    
      } else {
        # use character
        .JavaMethod("SJDBCBridge", "sjdbcAddStringColumn", "(Ljava/lang/String;Ljava/lang/String;[Ljava/lang/String;)V",
                    dfid, names(data)[i], as.character(data[,i]))
      }
    }
  }

  export.start <- function(dfid, driverClass, con, user, password,
                           table, appendToTable, preserveColumnCase) {
    .JavaMethod("SJDBCBridge", "sjdbcExportStart", 
                "(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;ZZ)V",
                dfid, driverClass, con, user, password, table, appendToTable, preserveColumnCase)
  }

  export.data <- function(dfid, commit=FALSE, batchSize) {
    .JavaMethod("SJDBCBridge", "sjdbcExportData", "(Ljava/lang/String;ZI)I", dfid, commit, batchSize)
  }

  # clear saved data
  on.exit(send.data.to.sjdbc(dfid, NULL))
  

	send.data.to.sjdbc(dfid, data)
	export.start(dfid, driverClass, con, user, password, table, appendToTable, preserveColumnCase)
	res <- export.data(dfid, TRUE, batchSize)

	if (!keepAlive) sjdbcCloseConnection()
	res
}
