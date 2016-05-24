epOutputHandler <-
function(res=NULL,epPlotInfo=NULL){
	
	if(!is.null(res) && !is.null(epPlotInfo)){
		final.output <- list(ExPosition.Data=res,Plotting.Data=epPlotInfo)
		class(final.output) <- c("expoOutput","list")	
		return(final.output)
	}#else if(!is.null(res) && is.null(epPlotInfo)){
	#	return(res)
	#}
	else{
		print("Unknown inputs. epOutputHandler must exit.")	
		return(0)
	}
	
	print("It is unknown how this was executed. epOutputHandler must exit.")
	return(0)
}
