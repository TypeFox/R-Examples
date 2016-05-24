tepOutputHandler <-
function(res=NULL,tepPlotInfo=NULL){
	
	if(!is.null(res) && !is.null(tepPlotInfo)){
		final.output <- list(TExPosition.Data=res,Plotting.Data=tepPlotInfo)
		class(final.output) <- c("texpoOutput","list")	
		return(final.output)
	}else if(!is.null(res) && is.null(tepPlotInfo)){
		return(res)
	}else{
		print("Unknown inputs. tepOutputHandler must exit.")	
		return(0)
	}
	
	print("It is unknown how this was executed. tepOutputHandler must exit.")
	return(0)
}
