distributionplot.POSIXct <- function(values, ...){
	# remove date
	# values <- strptime(format(values, "%H:%M:%S"), format="%H:%M:%S");
	myplot <- qplot(values, ...);
	return(myplot);	
}

distributionplot.hours_before_now <- function(values, ...){
	myfactor <- as.vector(values);
	distributionplot.do(myfactor, ...);	
}

distributionplot.multifactor <- function(values, ...){
	myfactor <- as.vector(values);
	distributionplot.do(myfactor, ...);
}

distributionplot.numeric <- function(values, ...){
	
	#exception if there are only a couple of unique values:
	if(length(unique(values)) < 8){
		values <- factor(values, ordered=T);
		return(distributionplot.do(values, ...));
	}
	
	myplot <- qplot(values, geom="bar", ...) 
	return(myplot);
}

distributionplot.factor <- function(values, ...){
	myplot <- qplot(values, geom="bar", fill=values, ...);
	if(length(levels(values)) > 7){
		myplot <- myplot + opts(axis.text.x=theme_text(angle=45));
	}	
	return(myplot);
}

distributionplot.character <- function(values, ...){
	#some string manipulation
	bigstring <- paste(values, collapse=" ");
	bigstring <- gsub("[\f\t.,;:`'\\\"\\(\\)<>]+", " ", bigstring);
	allwords <- tolower(strsplit(bigstring, " +")[[1]])
	
	#count, sort, head
	#library(reshape); #reshape::melt does not work.
	cloud <- melt(head(as.list(sort(table(allwords), decreasing=TRUE)), 200));
	words <- cloud[[2]];
	freq <- cloud[[1]];
	
	#make plot
	#myplot <- qplot(x=runif(length(words)), y=runif(length(words)), ..., geom="text", label=words, size=freq, color=freq) +
	#scale_size(range = c(6, 12)) + opts(axis.text.x = theme_blank()) + opts(axis.text.y = theme_blank()); 
	#return(myplot);
	
	#make plot with wordcloud package
	pal <- c("#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B")
	suppressWarnings(wordcloud(words, ceiling(log(freq+1,2)), min.freq=1, c(5,.5), random.order=FALSE, colors=pal))
	invisible()
	
	#myplot <- qplot(x=runif(length(words)), y=runif(length(words)), ..., geom="text", label=words, size=freq, color=freq) +
	#scale_size(range = c(6, 12))
	#return(myplot);
	
}

distributionplot.do <- function(values, ...){
	UseMethod("distributionplot");
}

#' Shows a histogram or barchart of the data
#' @param campaign_urn campaign id 
#' @param prompt_id id of the prompt
#' @param ... other arguments passed to oh.survey_response.read
#' @return ggplot2 plot object
#' @import wordcloud
#' @export
distributionplot <- function(campaign_urn, prompt_id, ...){
	
	#secret argument printurl for debugging	
	geturl(match.call(expand.dots=T));
		
	#get data
	myData <- oh.survey_response.read(campaign_urn=campaign_urn, prompt_id_list=prompt_id, column_list="urn:ohmage:prompt:response", ...);
	if(nrow(myData) > 0) myData <- na.omit(myData);
	fullname <- paste("prompt.id.", prompt_id, sep="");

	#check for no data
	if(nrow(myData) == 0 || sum(!is.na(myData[[fullname]])) == 0){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab="", ylab=""));
	}	
	
	#draw plot
	plottitle <- paste("distributionplot: ", prompt_id, sep="");	
	myplot <- distributionplot.do(na.omit(myData[[fullname]]), xlab="", ylab="", main=plottitle);
	 
	#return
	return(myplot);
}