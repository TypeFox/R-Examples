timeplot.POSIXct <- function(values, dates, aggregate, ...){
	#remove date
	newvalues <- as.numeric(format(values, "%H"));
	timeplot.do(newvalues, dates=dates, aggregate=aggregate, ...);
}

timeplot.hours_before_now <- function(values, dates, aggregate, ...){
	#same plot as timestamp but with hours subtracted 
	newvalues <- (as.numeric(format(dates, "%H")) - values + 24) %% 24;
	timeplot.do(newvalues, dates=dates, aggregate=aggregate, ...);
}


timeplot.multifactor <- function(values, dates, aggregate, ...){
	newvalues <- as.vector(values);
	newdates <- rep(dates, facdim(values));
	timeplot.do(newvalues, newdates, aggregate,...);
}

timeplot.numeric <- function(values, dates, aggregate, ...){

	dates <- as.Date(dates);
	mybinwidth <- aggregate;
	
	bindata <- bin.by.date(dates, values, binwidth=mybinwidth, probs=c(0,.5,1));
	myplot <- qplot(x=Date, y=Mean, ymin=Min, ymax=Max, data=bindata, ...) +
	geom_ribbon(alpha=0.5) +
	geom_line(size=1, color="blue") +
	geom_point(size=3, color="red");
     
	return(myplot);
}

timeplot.factor <- function(values, dates, aggregate, ...){
	
	dates <- as.Date(dates);
	mybinwidth <- aggregate;
	
	myData <- data.frame(date=dates, value=values);
	myplot <- ggplot(aes(x=date, fill=value), data=myData) + geom_bar(binwidth=mybinwidth, color="white", size=0.2);
	return(myplot);		
}

timeplot.character <- function(values, dates, ...){
	
	#create dataframe of strings
	dates <- as.Date(dates);
	#dates <- factor(unclass(dates), levels=seq(min(dates), max(dates), by=1));
	y <- runif(length(dates),0,1);
	angle <- rnorm(length(dates),0,10)
	myData <- data.frame(date=dates, text=values, y=y, angle=angle);
	
	#create plot
	myplot <- qplot(date, y, label=text, angle=angle, geom="text", data=myData, ...)
	return(myplot);
	
}

timeplot.default <- function(values, dates, ...){
	stop("No timeplot has been defined for variables of class: ", class(values))
}


timeplot.do <- function(values, dates, ...){
	UseMethod("timeplot")	
}

timeplot_with_aggregate <- function(values, dates, aggregate, ...){
	dates <- as.Date(dates);
	if(missing(aggregate) || is.null(aggregate)){
		totalperiod <- unclass(range(dates)[2] - range(dates)[1]);
		if(totalperiod < 30){
			aggregate <- 1;
		} else if (totalperiod < 180 ){
			aggregate <- 7;
		} else {
			aggregate <- 30;
		}
	} else {
		if(!is.numeric(aggregate)){
			stop("Argument aggregate has to be a number that represents the number of days to aggregate over.")
		}
	}
	
	timeplot.do(values, dates, aggregate, ...);
}

#' Timeseries plot of a prompt 
#' @param campaign_urn campaign id
#' @param prompt_id prompt id
#' @param aggregate number of days to aggregate over. Defaults to something smart.
#' @param ... other arguments passed on to oh.survey_response.read
#' @return ggplot2 plot object
#' @export
timeplot <- function(campaign_urn, prompt_id, aggregate, ...){
	
	#printurl
	geturl(match.call(expand.dots=T));
	
	#get data
	myData <- oh.survey_response.read(campaign_urn, prompt_id_list=prompt_id, column_list="urn:ohmage:prompt:response,urn:ohmage:context:timestamp", ...);
	if(nrow(myData) > 0) myData <- na.omit(myData);
	fullname <- paste("prompt.id.", prompt_id, sep="");
	
	#check for no data
	if(nrow(myData) == 0 || all(is.na(myData[[fullname]]))){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab="", ylab=""));
	}	

	#set aggregate
	if(missing(aggregate)){
		aggregate <- NULL;
	}
	
	#draw plot
	plottitle <- paste("timeplot: ", prompt_id, sep="");	
	myplot <- timeplot_with_aggregate(myData[[fullname]], myData$context.timestamp, aggregate=aggregate) + opts(title=plottitle) + xlab("") + ylab("");
	
	#return
	return(myplot);
}

