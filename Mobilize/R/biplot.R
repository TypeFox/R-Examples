biplot.POSIXct <- function(xvar, yvar, ...){
	if(is.character(yvar)){
		return(biplot.do(yvar, xvar, ...) + coord_flip());
	}
	myplot <- timeplot.do(yvar, xvar, ...);
	return(myplot);
}

biplot.character <- function(xvar, yvar, ...){
	if(is.character(yvar)){
		myplot <- distributionplot.do(c(xvar,yvar), ...);
	} else if("POSIXt" %in% class(yvar)){
		myplot <- timeplot.do(xvar, yvar, ...) + coord_flip();
	} else {
		x <- runif(n=length(yvar), 0, 1);
		myplot <- qplot(x=x, y=yvar, ...) + geom_text(aes(label=xvar)) + opts(axis.text.x = theme_blank()); 
	}
}

biplot.numeric <- function(xvar, yvar, xlab, ylab, ...){
	if("factor" %in% class(yvar)){
		#we switch the labels and axes
		myplot <- qplot(yvar, xvar, geom="boxplot", xlab=ylab, ylab=xlab, ...) + coord_flip();
		return(myplot);
	} else if(is.character(yvar)){
		return(biplot.do(yvar, xvar, ...) + coord_flip());
	} else {
		myplot <- qplot(xvar, yvar, geom="point", xlab=xlab, ylab=ylab, ...) + stat_density2d(aes(color = ..level..), geom="density2d")
		return(myplot);
	}
}

biplot.factor <- function(xvar, yvar, ...){
	if("factor" %in% class(yvar)){
		return(biplot.factorfactor(xvar,yvar, ...));
	} else if(is.character(yvar)){
		return(biplot.do(yvar, xvar, ...) + coord_flip());
	} else {
		myplot <- qplot(xvar, yvar, geom="boxplot", ...)
		return(myplot);
	}
}

biplot.factorfactor <- function(xvar, yvar, ...){
	#melt data into df
	#library(reshape);
	myData <- melt(table(xvar,yvar));
	names(myData) <- c("xvar", "yvar", "value");
	myData$xvar <- factor(myData$xvar, levels=levels(xvar), ordered=T);
	myData$yvar <- factor(myData$yvar, levels=levels(yvar), ordered=T);
	myData <- myData[myData$value > 0,];
	
	#make plot
	myplot <- qplot(x=xvar, y=yvar, size=value*2, color=value, label=value, data=myData, ...) + geom_point() +
	geom_text(aes(size=value), color="white") +
	scale_size(range = c(5, 20), guide="none");	

	#return plot 
	return(myplot);
}

biplot.do <- function(values, dates, ...){
	UseMethod("biplot")	
}


#' Generate a biplot of two variables
#' @param campaign_urn id of the campaign
#' @param prompt_id prompt on the x axis
#' @param prompt2_id prompt on the y axis
#' @param ... other parameters passed on to oh.survey_response/read
#' @return ggplot2 plot object
#' @export
biplot <- function(campaign_urn, prompt_id, prompt2_id, ...){
	
	#secret argument printurl for debugging
	geturl(match.call(expand.dots=T));
	
	#get data for both prompts
	myData <- oh.survey_response.read(campaign_urn, column_list="urn:ohmage:prompt:response", prompt_id_list=unique(c(prompt_id, prompt2_id)), ...);
	myData <- na.omit(myData);
	
	#check for empty plot
	if(nrow(myData) == 0){
		return(qplot(0,0,geom="text", label="request returned no data.", xlab=prompt_id, ylab=prompt2_id));
	}	
	
	xvarname <- paste("prompt.id.", prompt_id, sep="");
	yvarname <- paste("prompt.id.", prompt2_id, sep="");	
	
	x_var <- myData[[xvarname]];
	y_var <- myData[[yvarname]];
	
	#expand multifactors
	if(is.multifactor(x_var)){
		y_var <- rep(y_var, facdim(x_var));
		x_var <- as.vector(x_var);
	}
	
	if(is.multifactor(y_var)){
		x_var <- rep(x_var, facdim(y_var));
		y_var <- as.vector(y_var);
	}		
	
	#check for characters (multifactors have been converted already)
	#if(is.character(x_var)){
	#	stop("Biplots do not support text prompts for now. (", xvarname, ")\n")
	#}
	#
	#if(is.character(y_var)){
	#	stop("Biplots do not support text prompts for now. (", yvarname, ")\n")
	#}	
	
	myplot <- biplot.do(x_var, y_var, xlab=prompt_id, ylab=prompt2_id, main="");
	
	#return
	return(myplot);

}

