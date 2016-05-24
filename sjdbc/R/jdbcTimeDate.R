"jdbcTimeDate" <- 
function(data) {
	# converts a vector of dates data into JDBC standard Timestamp format
	# yyyy-mm-dd hh:mm:ss.fffffffff    (in GMT)
	format <- "%Y-%m-%d %H:%M:%OS3"
	#data@time.zone <- "GMT"
	as.character(data, format=format)
}