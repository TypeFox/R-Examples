responseplot.do <- function(dates, surveyvec, aggregate, ...){
	dates <- as.Date(dates);
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
	
	myData <- data.frame(date=dates, survey=surveyvec);
	myplot <- ggplot(aes(x=date, fill=survey), data=myData) + geom_bar(binwidth=mybinwidth);
	return(myplot);	
}



#' Create a responseplot
#' @param campaign_urn id of the campaign
#' @param aggregate optional number of days to aggregate over. Defaults to something smart.
#' @param privacy_state either "shared" or "private" or "both"
#' @param ... stuff to pass on to oh.survey_response.read
#' @return a responseplot
#' @import Ohmage ggplot2 methods
#' @export
#' @examples library(Mobilize)
#' \dontrun{
#' #authentication works like a cookie.
#' #oh.login("ohmage.admin", "ohmage.passwd", "https://example.com/app")
#' 
#' #list campaigns you are in
#' #oh.campaign.read()
#' 
#' #make a plot
#' #responseplot("urn:ohmage:campaign:mycampaign");
#' }
responseplot <- function(campaign_urn, aggregate, privacy_state="both", ...){
	
	#printurl
	geturl(match.call(expand.dots=T));
	
	#grab data
	myData <- oh.survey_response.function.read(campaign_urn, ...);
	if(nrow(myData) > 0) myData <- na.omit(myData);

	#filter shared.
	if(privacy_state == "shared"){
		myData <- myData[myData$privacy_state == "shared",];
	} else if (privacy_state == "private"){
		myData <- myData[myData$privacy_state == "private",];
	}	
	
	#check for no data
	if(nrow(myData) == 0){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab="", ylab=""));
	}	

	#draw plot
	plottitle <- paste("responseplot: ", gsub("urn:campaign:","",campaign_urn), sep="");
	myplot <- responseplot.do(
			rep(myData$date, myData$count), rep(myData$survey_id, myData$count),  
			xlab="", ylab="Response Count", main=plottitle, aggregate, ...);
	
	#return
	return(myplot)	
}

