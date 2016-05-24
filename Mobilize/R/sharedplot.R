#' A barchart of the number of shared and unshared responses per campaign
#' @param campaign_urn campaign id
#' @param ... arguments passed on to oh.survey_response.read
#' @return a ggplot2 object 
#' @export
sharedplot <- function(campaign_urn, ...){
	
	#check for secret 'printurl' argument
	geturl(match.call(expand.dots=T));
	
	#get data
	myData <- oh.survey_response.function.read(campaign_urn, privacy_state_item_list="survey", ...);
	
	#check if we have some data
	if(nrow(myData) == 0){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab="", ylab=""));
	}
		
	#make plot
	plottitle <- paste("sharedplot: ", gsub("urn:campaign:","",campaign_urn), sep="");
	
	myplot <- ggplot(aes(x=survey_id, y=count, group=privacy_state, fill=privacy_state), data=myData) +
		geom_bar(stat="identity") # (xlab="Survey", ylab="count", main=plottitle);

	#return
	return(myplot);
}