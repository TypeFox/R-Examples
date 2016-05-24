fun.classCheck <-
function(YPmodelResult,...)
{
	
	if(class(YPmodelResult) != 'YPmodel'){
		stop("Error: ", paste(sQuote("YPmodelResult"), sep = ""), " must be of class ", 
	            paste(dQuote("YPmodel"), sep = ""), sep = "")
	}

}
