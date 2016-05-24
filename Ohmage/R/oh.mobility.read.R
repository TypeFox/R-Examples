#' Read Mobility data
#' 
#' @param date date in ISO format
#' @param username name of user to query (only works when server allows to see others data). 
#' @param column_list variables to be returned
#' @param ... stuff passed to oh.call
#' @return a dataframe with mobility data
#' @export
oh.mobility.read <- function(date = today(), username=getOption("ohmage_username"), column_list, ...){
	if(missing(column_list)){
		column_list <- c("mobility:id", "location", "mobility:mode", "mobility:timestamp");
	}
	if(is.character(date) && nchar(date) != "10"){
		stop("Date has to be in format YYYY-mm-dd");
	} 
	if("Date" %in% class(date)){
		date <- as.character(date);
	}
	if("POSIXt" %in% class(date)){
		date <- as.character(as.Date(date));
	}

	xhr <- oh.call("/mobility/read", date=date, username=username, column_list=paste(column_list, collapse=","), ...);		
	
	output <- as.data.frame(do.call("rbind",lapply(lapply(lapply(xhr$data, "[[", "l"), parsevector), unlist)), stringsAsFactors=FALSE);
	if(nrow(output) == 0){
		output <- as.data.frame(matrix(nrow=length(xhr$data), ncol=0));
	}
	output$id <- unlist(lapply(xhr$data, "[[", "id"));
	output$m  <- unlist(lapply(xhr$data, "[[", "m"));
	output$ts <- unlist(lapply(xhr$data, "[[", "ts"));
	output$lo <- as.numeric(output$lo);
	output$la <- as.numeric(output$la);
	output$ac <- as.numeric(output$ac);
	output$ts <- strptime(output$ts, format="%Y-%m-%d %H:%M:%S");

	if("mobility:sensor_data" %in% column_list){
		output$t  <- structure(as.numeric(unlist(lapply(xhr$data, "[[", "t")))/1000, class=class(Sys.time()));		
		output$speed <- as.numeric(unlist(lapply(lapply(xhr$data, "[[", "data"), "[[", "sp")));
	}
	
	#sort
	output <- output[order(output$ts),];
	return(output);
}
