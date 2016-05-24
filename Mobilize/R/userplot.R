userplot.POSIXct <- function(values, dates, ...){
	myplot <- timeplot_with_aggregate(values, dates, ...);
	return(myplot);
}

userplot.hours_before_now <- function(values, dates, ...){
	myplot <- timeplot_with_aggregate(values, dates, ...);
	return(myplot);
}

userplot.multifactor <- function(values, dates, ...){
	newvalues <- as.vector(values);
	newdates <- rep(dates, facdim(values));
	userplot.do(newvalues, newdates, ...);
}

userplot.numeric <- function(values, dates, ...){
	myplot <- timeplot_with_aggregate(values, dates, ...);
	return(myplot);
}

userplot.factor <- function(values, dates, ...){
	dates <- as.Date(dates);
	dates <- factor(unclass(dates), levels=seq(min(dates), max(dates), by=1));

	#library(reshape);
	myData <- melt(table(dates, values));
	names(myData) <- c("dates", "values", "count");
	myData$dates <- as.Date(myData$dates);
	myData <- myData[myData$count > 0,];
	
	myplot <- qplot(x=dates,y=values, size=count*2, color=count, label=count, data=myData, ...) + geom_point() +
	geom_text(aes(size=count), color="white") +
	scale_size(range = c(5, 20), guide="none");

	return(myplot);
}

userplot.character <- function(values, dates, ...){
	
	#same as timeplot
	myplot <- timeplot_with_aggregate(values, dates, ...);
	return(myplot);
	
}

userplot.default <- function(values, dates, ...){
	stop("No userplot has been defined for variables of class: ", class(values))
}


userplot.do <- function(values, dates, ...){
	UseMethod("userplot")	
}

#' Timeseries plot of data for a single user
#' @param campaign_urn campaign id
#' @param prompt_id prompt id
#' @param user_id user id
#' @param ... arguments passed on to oh.survey_response.read
#' @return a ggplot2 plot object
#' @import reshape2
#' @export
userplot <- function(campaign_urn, prompt_id, user_id, ...){
	
	#printurl
	geturl(match.call(expand.dots=T));
	
	#get data
	myData <- oh.survey_response.read(campaign_urn, prompt_id_list=prompt_id, user_list=user_id, ...);
	myData <- na.omit(myData);
	fullname <- paste("prompt.id.", prompt_id, sep="");
	
	#check for now data
	if(nrow(myData) == 0 || all(is.na(myData[[fullname]]))){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab="", ylab=""));
	}	
	
	#make plot	
	plottitle <- paste("userplot: ", user_id, sep="");	
	myplot <- userplot.do(myData[[fullname]], myData$context.timestamp, xlab="", ylab=prompt_id, main=plottitle)
	
	#return
	return(myplot);
}


