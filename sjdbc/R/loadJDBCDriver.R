loadJDBCDriver <- function(driverJar) {
	.JavaMethod("SJDBCBridge", "sjdbcAddDrivers", "([Ljava/lang/String;)V", 
		driverJar)
	invisible()
}
