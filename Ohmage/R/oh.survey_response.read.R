#' Read survey responses
#' @param campaign_urn campaign id
#' @param prompt_id_list comma seperated ist of promt id's
#' @param privacy_state shared or private
#' @param column_list ohmage columns
#' @param output_format should be json-columns
#' @param user_list list of users
#' @param to.data.frame if data should be converted to a dataframe
#' @param ... other args passed to ohmage
#' @export
oh.survey_response.read <- function(campaign_urn, prompt_id_list="urn:ohmage:special:all", privacy_state="both", 
	column_list="urn:ohmage:user:id,urn:ohmage:prompt:response,urn:ohmage:context:timestamp", output_format="json-columns", 
	user_list="urn:ohmage:special:all", to.data.frame=TRUE, ...){

	#if prompt_id_list is a vector:
	prompt_id_list <- paste(prompt_id_list, collapse=",");
	column_list <- paste(column_list, collapse=",");
	user_list <- paste(user_list, collapse=",");
	
	if(privacy_state=="both"){
		xhr <- oh.call("/survey_response/read ",campaign_urn=campaign_urn, prompt_id_list=prompt_id_list, column_list=column_list, output_format=output_format, user_list=user_list, ...);			
	} else {	
		xhr <- oh.call("/survey_response/read ",campaign_urn=campaign_urn, prompt_id_list=prompt_id_list, column_list=column_list, output_format=output_format, user_list=user_list, privacy_state=privacy_state, ...);	
	}	

	if(!to.data.frame){
		return(xhr);
	}
	
	rows <- xhr$metadata$number_of_surveys;
	columns <- length(xhr$metadata$items);
	
	if(rows==0){
		#warning("This query returned no data.");
		itemnames <- unlist(xhr$metadata$items);
		emtpydf <-  as.data.frame(matrix(nrow=0, ncol=length(itemnames), dimnames=list(vector(),itemnames)));
		return(emtpydf);
	}
	
	if(length(xhr$data) != columns){
		stop("number of items in metadata does not match with data.")
	}
	
	mydf <- as.data.frame(matrix(NA, rows, 0));
	for(i in 1:columns){
		varname <- names(xhr$data[i]);
		varname <- gsub("urn:ohmage:", "", varname);
		varname <- gsub(":", ".", varname);
		
		#TODO: this check is sort of a hack for the no values bug in ohmage 2.3.
		newvar <- parse.item(xhr$data[i]);
		if(length(newvar) > 0){
			if(length(newvar) != rows){
				stop("JSON decoding error. Number of values for variable: ", varname, " (n=", length(newvar), ") does not match with metadata.number_of_surveys (", rows, ")");
			}
			mydf[[varname]] <- newvar;
		} else {
			mydf[[varname]] <- NA;
		}
	}
	return(mydf);
}
