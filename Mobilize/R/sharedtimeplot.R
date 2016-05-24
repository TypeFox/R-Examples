sharedtimeplot.do <- function(dates, surveyvec, sharedvec, aggregate, ...){

	#remove time
	if(missing(aggregate)){
		totalperiod <- unclass(range(dates)[2] - range(dates)[1]);
		if(totalperiod < 30){
			mybinwidth <- 1;
		} else if (totalperiod < 180 ){
			mybinwidth <- 7;
		} else {
			mybinwidth <- 30;
		}
	} else {
		if(!is.numeric(aggregate)){
			stop("Argument aggregate has to be a number that represents the number of days to aggregate over.")
		}
		mybinwidth <- aggregate;
	}
	
	myData <- data.frame(date=dates, survey=surveyvec, privacy=sharedvec);
	myplot <- ggplot(aes(x=date, fill=privacy, group=privacy), data=myData) + geom_bar(binwidth=mybinwidth) + facet_wrap(~survey, ncol=1);
	return(myplot)
}



#' Timeseries plot of the number of shared and unshared responses per campaign
#' @param campaign_urn campaign id
#' @param aggregate number of days to aggregate over. Optional. Defaults to something smart.
#' @param ... other arguments passed to oh.survey_response.read
#' @return a ggplot2 object
#' @export
sharedtimeplot <- function(campaign_urn, aggregate, ...){
	
	#printurl
	geturl(match.call(expand.dots=T));
	
	#grab data
	myData <- oh.survey_response.function.read(campaign_urn, ...);
	if(nrow(myData) > 0) myData <- na.omit(myData);
	
	#check for no data
	if(nrow(myData) == 0){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab="", ylab=""));
	}	
	
	#draw plot
	plottitle <- paste("sharedtimeplot: ", gsub("urn:campaign:","",campaign_urn), sep="");
	myplot <- sharedtimeplot.do(rep(myData$date, myData$count), rep(myData$survey_id, myData$count), rep(myData$privacy_state, myData$count), aggregate, ...);

	return(myplot)
}
