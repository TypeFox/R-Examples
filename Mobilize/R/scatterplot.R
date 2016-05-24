#' Create a scatterplot of two prompts
#' @param campaign_urn id of the campaign
#' @param prompt_id prompt on the x axis.
#' @param prompt2_id prompt on the y axis
#' @param jitter T/F. If points should be jittered 
#' @param ... arguments passed on to oh.survey_response.read
#' @return ggplot2 object
#' @export
scatterplot <- function(campaign_urn, prompt_id, prompt2_id, jitter=TRUE, ...){
	
	#check for secret 'printurl' argument
	geturl(match.call(expand.dots=T));
	
	#get data
	myData <- oh.survey_response.read(campaign_urn, column_list="urn:ohmage:prompt:response", prompt_id_list=unique(c(prompt_id, prompt2_id)), ...);
	myData <- na.omit(myData);

	#check for empty plot
	if(nrow(myData) == 0){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab=prompt_id, ylab=prompt2_id));
	}		
	
	#expand multifactors
	xvarname <- paste("prompt.id.", prompt_id, sep="");
	yvarname <- paste("prompt.id.", prompt2_id, sep="");		
	x_var <- myData[[xvarname]];
	y_var <- myData[[yvarname]];
	
	if(is.multifactor(x_var)){
		y_var <- rep(y_var, facdim(x_var));
		x_var <- as.vector(x_var);
	}
	
	if(is.multifactor(y_var)){
		x_var <- rep(x_var, facdim(y_var));
		y_var <- as.vector(y_var);
	}		
	
	#draw plot
	myplot <- qplot(x_var, y_var, xlab=prompt_id, ylab=prompt2_id, main="", geom="blank");

	#if jitter == TRUE, add a little jitterish
	if(isTRUE(jitter)){
		myplot <- myplot + geom_jitter(position=position_jitter(width=.25, height=.25));		
	} else {
		myplot <- myplot + geom_point();		
	}
	
	if(is.factor(myData[[xvarname]]) && length(levels(myData[[xvarname]])) > 7){
		myplot <- myplot + opts(axis.text.x=theme_text(angle=45));
	}
	
	return(myplot);
}
