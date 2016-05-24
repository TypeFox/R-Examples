lhc_netlogo_graphMeasuresForParameterChange <- function(FILEPATH,PARAMETERS,MEASURES,MEASURE_SCALE,CORCOEFFSOUTPUTFILE,LHCSUMMARYFILENAME,TIMEPOINTS,TIMEPOINTSCALE)
{
	# Not using this anymore, as using check.names=FALSE when reading CSV files
	# Check the parameter and measures strings
	#PARAMETERS<-table_header_check(PARAMETERS)
	#MEASURES<-table_header_check(MEASURES)

	# Call the spartan function
	lhc_graphMeasuresForParameterChange(FILEPATH,PARAMETERS,
		MEASURES,MEASURE_SCALE,CORCOEFFSOUTPUTFILE,
		LHCSUMMARYFILENAME,TIMEPOINTS,TIMEPOINTSCALE)

}
