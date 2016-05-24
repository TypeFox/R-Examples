#' Survey Response Function Read
#' 
#' Allows users with the appropriate permissions to read aggregate data about survey responses. The API takes an id parameter that tells the system which function to perform on the data.
#' 
#' @param campaign_urn campaign urn
#' @param id Function id. The currently supported id is privacy_state
#' @param privacy_state depricated. Should be NULL. 
#' @param privacy_state_item_list comma seperated list of output
#' @param ... other stuff passed to oh.call
#' @return a list
#' @export
oh.survey_response.function.read <- function(campaign_urn, id="privacy_state", privacy_state=NULL, privacy_state_item_list="date,survey", ...){
	
	if(!is.null(privacy_state)){
		stop("Privacy_state argument was matched for function read. This should never happen. Probably conflicting with privacy_state_item_list.")
	}
	
	xhr <- oh.call("/survey_response/function/read", campaign_urn=campaign_urn, id=id, privacy_state_item_list=privacy_state_item_list, ...)
	mydata <- xhr$data;
	
	
	if(length(mydata) == 0){
		return(data.frame());
	}
	
	output <- list();
	for(i in 1:length(mydata)){
		output[[i]] <- cbind(do.call(rbind,lapply(mydata[[i]], unlist)), idfunname=names(mydata[i]));
	}
	output <- as.data.frame(do.call(rbind, output), stringsAsFactors=FALSE);
	names(output)[which(names(output)=="idfunname")] <- id;
	
	if(!is.null(output$date)){
		output$date <- as.Date(output$date);
	}

	output$count <- as.numeric(output$count);
	return(output);
}
