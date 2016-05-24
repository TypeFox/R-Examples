oat_generate_netlogo_behaviour_space_XML <-
function(FILEPATH,NETLOGO_SETUPFILE_NAME,PARAMETERS,PARAMVALS,NETLOGO_SETUP_FUNCTION,NETLOGO_RUN_FUNCTION,MEASURES,EXPERIMENT_REPETITIONS,RUNMETRICS_EVERYSTEP)
{
	if(requireNamespace("XML",quietly=TRUE))
	{
		# START A NEW XML FILE, WITH EXPERIMENTS AS THE TOP TAG (AS REQUIRED BY NETLOGO)
		xml<-XML::xmlOutputDOM(tag="experiments")

		# NEXT TAG IN IS EXPERIMENT
		xml$addTag("experiment", attrs=c(name="OAT_Sample",repetitions=EXPERIMENT_REPETITIONS, runMetricsEveryStep=RUNMETRICS_EVERYSTEP),close=FALSE)

		# NOW THE PROCEDURES TO CALL SETUP, GO, AND OUTPUT MEASURES TO ANALYSE
		xml$addTag("setup",NETLOGO_SETUP_FUNCTION)
		xml$addTag("go",NETLOGO_RUN_FUNCTION)

		## NOW TO DO THE MEASURES
		for(MEASURE in 1:length(MEASURES))
		{
			xml$addTag("metric",MEASURES[MEASURE])	
		}

		for(PARAM in 1:length(PARAMETERS))
		{
			# NOW SOME PARAMETERS ARE BEING VARIED, SOME NOT
			# THE ONES THAT ARE BEING VARIED HAVE A MIN, MAX, AND INCREMENT IN SQUARE BRACKETS, SEPARATED BY A COMMA
			# THUS WE CAN DISTINGUISH THESE WITH A STRING TOKENIZER
			PARAMVALSPLIT<-(strsplit(PARAMVALS[PARAM],","))[[1]]

			if(length(PARAMVALSPLIT)==1)
			{
				# THIS PARAMETER IS NOT BEING VARIED, AND THUS WE CAN JUST SPECIFY THE PARAMETER VALUE
				# FOR EACH PARAMETER, ADD THE ENUMERATEDVALUESET TAG, SIMULATION VARIABLE NAME
				xml$addTag("enumeratedValueSet", attrs=c(variable=(PARAMETERS[PARAM])),close=FALSE)
	
				# NOW ADD THE VALUE
				xml$addTag("value", attrs=c(value=(PARAMVALS[PARAM])))

				# CLOSE THE ENUMERATED VALUE SET TAG
				xml$closeTag()
			}else
			{
				# THIS IS A PARAMETER BEING ANALYSED, AND THUS WE NEED TO TELL NETLOGO TO ALTER THE VALUES
				# GET THE MIN, MAX, AND INCREMENT
				# NOTE FOR MIN AND INCREMENT, WE NEED TO REMOVE THE OPENING AND CLOSING SQUARE BRACKET
				MIN<-substring(PARAMVALSPLIT[[1]],2)
				MAX<-PARAMVALSPLIT[[2]]
				INC<-substring(PARAMVALSPLIT[[3]], 1, nchar(PARAMVALSPLIT[[3]])-1)

				# NOW TO BUILD THE TAGS
				xml$addTag("steppedValueSet", attrs=c(variable=(PARAMETERS[PARAM]),first=MIN,step=INC,last=MAX))
			}
		}


		# CLOSE THE EXPERIMENT TAG
		xml$closeTag()	

		# CLOSE THE EXPERIMENTS TAG
		xml$closeTag()	

		XML::saveXML(xml,file=paste(FILEPATH,"/",NETLOGO_SETUPFILE_NAME,".xml",sep=""),indent=TRUE, prefix = '<?xml version="1.0" encoding="us-ascii"?>\n',
		doctype = '<!DOCTYPE experiments SYSTEM "behaviorspace.dtd">')
	}
	else
	{
		print("The oat_generate_netlogo_behaviour_space_XML function requires the XML package to be installed")
	}
}
