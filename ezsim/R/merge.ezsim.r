#' Merge two ezsim objects. Either \code{m} or \code{parameter_def} of two ezsim objects must be the same. If \code{parameter_def} are the same, the merging is regarded as increasing the number of simulation \code{m}. If the \code{parameter_def} are different and \code{m} are the same, the merging is regarded as extending the \code{parameter_def}. Notice that, use_seed will not be merged!
#' @name merge.ezsim
#' @aliases merge.ezsim
#' @title Merge two ezsim objects
#' @method merge ezsim
#' @param x A ezsim to merge with 
#' @param y A ezsim to merge with
#' @param ... unused
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @seealso \code{\link{ezsim}}

merge.ezsim<-function(x,y,...){		
	# Check whether two ezsim objects are the same
	check_cond<-
	foreach (i = c('estimator','true_value','dgp','display_name'),.combine=c,.final=all) %do% {
		digest(x[[i]])==digest(y[[i]])
	}
	if (!check_cond)
		stop('Two ezsim objects are not the same')

	## Same parameter_def, adding the number of simulation. 
	if (digest(x$parameter_def)==digest(y$parameter_def)){
	
		out<-x
		out$m<-x$m+y$m
		
		i=j=NULL
		## merge the simulation_result
		out$simulation_result<-
		foreach (i = x$simulation_result, j = y$simulation_result ) %do% {
			c(i,j)
		}
		## reproduce the simulation table
		tryCatch({
			out<-createSimulationTable(out)
		}, error = function(e){
			cat("Error in creating simulation table\n")
			stop(e)
		})
		return(out)
	}
	## Same m, different parameter
	if ( (x$m == y$m) & (digest(x$parameter_def)!=digest(y$parameter_def)) ) {
		out<-x
		
		# check whether there's error when merging x$parameter_def and y$parameter_def	
		tryCatch({
			out$parameter_def<-merge(x$parameter_def,y$parameter_def)
		},error = function(e) {
			cat("Error in merging parameterDef")
			stop(e)
		})

		## merge sim
		warning("Merging simulation with parameterDef is not recommended. Potential problems : number of simulation for some cases (same set of parameter appear in both ezsim object) is doubled; Some cases is missing (The cardinality of outer product of two parameterDef is larger than concatenating of them). If you don't know what is going on, run the ezsim object again. e.g. run(ezim_object)")
		out$parameter_list<-c(x$parameter_list,y$parameter_list)
		out$simulation_result=c(x$simulation_result,y$simulation_result)
		out$simulation_table=rbind(x$simulation_table,y$simulation_table)
		
		return(out)		
	}
	stop("Conditions for merging is not satisfied. Please check")
}

# merge.ezsim<-function(x,y){
	# if (class(y)!='ezsim')
		# stop('y must be an ezsim object')
		
	# check_cond<-
	# foreach (i = c('estimator','true_value','dgp','display_name','packages'),.combine=c,.final=all) %do% {
		# digest(x[[i]])==digest(y[[i]])
	# }
	# if (!check_cond)
		# stop('Two ezsim objects are not the same')

	# # Same parameter_def, adding the number of simulation. 
	# if (digest(x$parameter_def)==digest(y$parameter_def)){
		# if (digest(x$TV_table)!=digest(y$TV_table))
			# stop('True value tables are not the same')
	
		# out<-x
		# out$m<-x$m+y$m
		# out$sim<-rbind(x$sim,y$sim)
		# return(out)
	# }
	# # Same m, different parameter
	# if ( (x$m == y$m) & (digest(x$parameter_def)!=digest(y$parameter_def)) ) {
	
		# out<-x
		# # check whether there's error when merging x$parameter_def and y$parameter_def	
		# out$parameter_def<-
		# tryCatch( merge( x$parameter_def,y$parameter_def) , error = function(e) stop(e) )
				
		# # merge sim
		# out$sim<-rbind(x$sim,y$sim)
		
		# # merge TV_table (regenerate it )
		# out$TV_table<-rbind(x$TV_table,y$TV_table)
		
		# estimator_name<-setdiff(names(out$TV_table),c('estimator','value'))
		
        # estimator_name<-estimator_name[length(estimator_name):1]
      
		# TV_table_order<-
		# eval(parse(text=paste('with(out$TV_table, order(','estimator,', paste(estimator_name,collapse=','),'))',sep='')))

		# out$TV_table<-out$TV_table[TV_table_order,]
				
		# return(out)		
	# }
# }
