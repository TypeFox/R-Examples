fun.errorMessage <-
function(type)
{
	if(type=='DataSet'){
		txt <- "The data sets was not specificated."
	}

	if(type=='DataLength'){
		txt <- "The length of data must be set."
	}


	if(type=='DefaultParameter'){
		txt <- "Parameters are not specificated, using the default one."
	}

	return(txt)

}
