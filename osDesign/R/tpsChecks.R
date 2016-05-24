tpsChecks <-
function(X,
										  strata,
											nII,
											cohort,
											NI,
											nII0=NULL,
											nII1=NULL,
											nCC=NULL,
											threshold=NULL)
{
  ##
  if(!is.list(strata))
  {
  	##
  	if(sum(!is.element(strata, 0:ncol(X))) > 0)
  		return("* 'strata' is invalid")
		##
  	if(max(strata) == 0)
  	{
  		if(is.null(nII))
  			return("* 'nII' is required when strata == 0")
			if(!is.null(nII0))
  			print("* Warning: argument 'nII0' is ignored when strata == 0")
			if(!is.null(nII1))
  			print("* Warning: argument 'nII1' is ignored when strata == 0")
  	}
	  ##
  	if(max(strata) > 0)
  	{
  		if(is.null(nII) & (is.null(nII0) | is.null(nII1)))
  			return("* Require valid phase II sample sizes: (i) 'nII' or (ii) 'nII0' and 'nII1'")
			if(is.element(0, strata))
  			print("* Warning: ignoring strata == 0")
  	}
  }

	##
  if(!is.null(nII))
  {
		if(min(nII) < 0 & length(nII) == 1)
			return("* Phase II sample size 'nII' is negative")
		if(min(nII) < 0 & length(nII) > 1)
			return("* At least one phase II sample size 'nII' is negative")
  }
  
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
  if(!is.null(nCC))
  {
		if(min(nCC) < 0 & length(nCC) == 1)
			return("* Case-control sample size 'nCC' is negative")
		if(min(nCC) < 0 & length(nCC) > 1)
			return("* At least one case-control sample size 'nCC' is negative")
  }

 	##
	if(!is.null(threshold))
	{
		if(length(threshold) != 2)
			return("* 'threshold' is not a pair of numbers")
	}
	
	##
	return("")
}
