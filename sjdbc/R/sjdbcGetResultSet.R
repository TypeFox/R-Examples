"sjdbcGetResultSet"<-
function(key, unregister = TRUE, default.num.rows = NULL,
         start.at.first=TRUE, rows.to.read=-1)
{
	# [Adapted from javaGetResultSet() in S-PLUS 8.0.4]
	#
	# Get the data from the ResultSet into the statics in 
	# SJDBCResultSetUtilities.  Returns NULL if key not recognized.
	#
	# key: String giving the key into the Hashtable in SJDBCResultSetUtilities.  The
	#      ResultSet needs to be put into this Hashtable in Java.
	# unregister: Should the ResultSet be removed from the Hashtable after the 
	#             data is retrieved.
	# default.num.rows: When the ResultSet is of type ResultSet.TYPE_FORWARD_ONLY, 
	#             the number of rows cannot be determined until after all of
	#             the data has been retrieved.  If the ResultSet has more than
	#             the default number, the array sizes will be doubled whenever the
	#             current capacity is reached.  If the ResultSet is not of
	#             TYPE_FORWARD_ONLY, this is not used.
	# start.at.first : if true, set the ResultSet to start with the first row before
	#             reading.  Otherwise, start with the current row.
	# rows.to.read : if less than zero, read all rows in the result set.  Otherwise,
	#             this is the maximum number of rows to read.
	className <- "SJDBCResultSetUtilities"
	if(!is.null(default.num.rows)) {
		.JavaMethod(className, "setDefaultNumberOfRows", "(I)V", 
			default.num.rows)
	}
	
	.JavaMethod(className, "snextPopulateData", 
		"(Ljava/lang/String;ZI)V", key, start.at.first, rows.to.read)
		
	colNames <- .JavaMethod(className, "snextGetColNames", "(Ljava/lang/String;)[Ljava/lang/String;", key) 
	colTypeCodes <- .JavaMethod(className, "snextGetColTypeCodes", "(Ljava/lang/String;)[I", key) 
	data <- data.frame()
	for (i in 1:length(colTypeCodes)) {
		if (colTypeCodes[i] %in% c(2,8)) 
			colData <- .JavaMethod(className, "snextGetDoubleColumn", "(Ljava/lang/String;I)[D", key, i-1)
		#else if (colTypeCodes[i] == 6) 
		#	colData <- .JavaMethod(className, "snextGetFloatColumn", "(Ljava/lang/String;I)[F", key, i-1)
		else if (colTypeCodes[i] == 4)
			colData <- .JavaMethod(className, "snextGetIntegerColumn", "(Ljava/lang/String;I)[I", key, i-1)
		else if (colTypeCodes[i] == -7) 
			colData <- .JavaMethod(className, "snextGetBooleanColumn", "(Ljava/lang/String;I)[Z", key, i-1)
		else  colData <- .JavaMethod(className, "snextGetStringColumn", "(Ljava/lang/String;I)[Ljava/lang/String;", key, i-1)
		
		if (nrow(data) == 0) data <- data.frame(colData) else data <- data.frame(data, colData)
	}
	
	names(data) <- colNames
	
	# Release the data
	.JavaMethod(className, "snextReleaseData", "(Ljava/lang/String;)V", key)
	
	# Release the ResultSet reference
	if(unregister) {
		.JavaMethod(className, "unregister", "(Ljava/lang/String;)V",
			key)
	}
	
	# Can use the colTypes to convert Date, Time, Timestamp, and
	# possibly other classes to similar S-PLUS classes.  
	# Uses strings to pass data back in standard JDBC format
	# due to loss-of-precision when passing back as doubles and
	# reconverting to julian/millisecond values.  
	## 91 - Date, 92 - Time, 93 - Timestamp
	pos <- match(colTypeCodes, c(91, 92, 93))
	pos <- (1:length(pos))[!is.na(pos)]
	if(length(pos) > 0) {
		for(i in 1:length(pos)) {
			data[, pos[i]] <- as.POSIXct(as.character(data[, pos[i]]), 
				format="%Y-%m-%d %H:%M:%S")
		}
	}
	return(data)
}
