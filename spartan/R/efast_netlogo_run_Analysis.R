efast_netlogo_run_Analysis<-function(FILEPATH,MEASURES,PARAMETERS,NUMCURVES,NUMSAMPLES,OUTPUTMEASURES_TO_TTEST,TTEST_CONF_INT,
						GRAPH_FLAG,EFASTRESULTFILENAME,TIMEPOINTS,TIMEPOINTSCALE)
{
	# Not using the header check anymore now we do not check the column names on read in.
	# Check the parameter and measures strings
	# PARAMETERS<-table_header_check(PARAMETERS)
	# MEASURES<-table_header_check(MEASURES)

	# Call the spartan function
	efast_run_Analysis(FILEPATH,MEASURES,PARAMETERS,NUMCURVES,NUMSAMPLES,OUTPUTMEASURES_TO_TTEST,TTEST_CONF_INT,
				GRAPH_FLAG,EFASTRESULTFILENAME,TIMEPOINTS,TIMEPOINTSCALE)
}
