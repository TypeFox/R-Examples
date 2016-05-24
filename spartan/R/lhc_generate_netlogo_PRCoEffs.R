lhc_generate_netlogo_PRCoEffs <- function(FILEPATH,PARAMETERS,MEASURES,LHCSUMMARYFILENAME,CORCOEFFSOUTPUTFILE)
{
	# Not using this anymore, as using check.names=FALSE when reading in CSV files
	# Check the parameter and measures strings
	#PARAMETERS<-table_header_check(PARAMETERS)
	#MEASURES<-table_header_check(MEASURES)

	# Call the spartan function
	lhc_generatePRCoEffs(FILEPATH,PARAMETERS,MEASURES,LHCSUMMARYFILENAME,CORCOEFFSOUTPUTFILE)
}
