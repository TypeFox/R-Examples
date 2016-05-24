#' Create a simulation table from simulation result 
#' @name createSimulationTable.ezsim
#' @aliases createSimulationTable.ezsim
#' @title createSimulationTable
#' @param x an ezsim object
#' @param \dots not used 
#' @export 
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com} 
#' @keywords internal
createSimulationTable.ezsim<-function(x,...){
	# for each parameter set, we have m simulation result, we will loop over them
	i=j=NULL
	x$simulation_table<-foreach(i = x$simulation_result , j = x$parameter_list,.combine=rbind ) %do% {
		
		# apply the estimator_parser
		out<-lapply(i,x$estimator_parser)
		
		if (any(!sapply(out,is.vector)))
			stop("estimator parser do not return a vector")
			
		## TODO built-in parser
		
		# if any of the estimator dont have a name, create a name for them.
		if (any(sapply(out,function(x) is.null(names(x)) ))){
			out<-
			foreach ( k = out ) %do% {
				if (is.null(names(k)))
					names(k)<-paste("estimates_",1:length(k))
				k
			}
		}
			
		## becareful about the name of estimator.
		## if ( is the starting value, may have problem
		## remove I() ?
		out<-data.frame(do.call(rbind,out))
		estimators_name <- names(out)
		out<-melt(out,measure.vars=estimators_name,variable_name='estimator')
		names(out)[2]<-"value_of_estimator"
		
		# compute the true value
		if (is.function(x$true_value)){
			true_value <-  data.frame(t(Jmisc::evalFunctionOnList(x$true_value,j)))
			names(true_value)<-estimators_name
			true_value <- melt(true_value,measure.vars=estimators_name,variable_name='estimator')
			names(true_value)[2]<-"value_of_TV"
			out <- merge(out,true_value,by="estimator",suffixes="")	
		} else {
			out$true_value=NA
		}
		
		par_value <- unlist(j[names(x$parameter_def$selection)])
		out <- Jmisc::addCol(out,value=par_value)
		out
	}
	return(x)
}
