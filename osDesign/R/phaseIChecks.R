phaseIChecks <-
function(X,
										     strata,
											   cohort,
											   NI,
											   nII0=NULL,
											   nII1=NULL)
{
  ##
  if(sum(!is.element(strata, 1:ncol(X))) > 0)
  	return("* 'strata' is invalid")

	##
 	if(cohort == FALSE)
 	{
 		if(is.null(NI))
  		return("* 'NI' must be specified if phase I arises via case-control sampling")
	  if(!is.null(NI))
  	{
  		if(length(NI) != 2)
  			return("* 'NI' should be a pair of Phase I sample sizes for controls and cases")
	  	if(min(NI) < 0)
  			return("* Phase I case-control sample size 'NI' is not positive")
  	}
 	}

	##
	return("")
}
