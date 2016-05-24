#' Score Profile Plot
#'
#' The \code{profileplot} function creates a profile plot for a matrix or dataframe with multiple scores or subscores using \code{\link[ggplot2]{ggplot}} function in \code{ggplot2} package.
#'
#' @export
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot aes element_blank element_rect geom_line geom_point scale_colour_hue scale_shape_discrete scale_x_discrete scale_y_continuous theme
#' @importFrom graphics plot points text
#' @importFrom RColorBrewer brewer.pal
#' @param form A matrix or dataframe including two or more subscores.
#' @param person.id A vector that includes person ID values (Optional).
#' @param standardize If not FALSE, all scores are rescaled with a mean of 0 and standard deviation of 1. Default is TRUE.
#' @param interval The number of equal intervals from the mimimum score to the meximum score. Default is 10. Ignored when by.pattern=FALSE.
#' @param by.pattern If TRUE, the function creates a profile plot with level and pattern values using ggplot2. Otherwise, the function creates a profile plot showing profile scores of persons using the base graphics in R. Default is TRUE.
#' @param original.names Use the original column names in the data. Otherwise, columns are renamed as v1,v2,.... Default is TRUE.
#' @return  The \code{profileplot} functions returns a score profile plot from either \link[ggplot2]{ggplot} or the base graphics in R.
#' @examples
#' \dontrun{
#'	data(PS)
#'  myplot <- profileplot(PS[,2:4], person.id = PS$Person,by.pattern = TRUE, original.names = TRUE)
#'  myplot
#' 
#' data(leisure)
#' leis.plot <- profileplot(leisure[,2:4],standardize=TRUE,by.pattern=FALSE)
#' leis.plot
#' }
#' @seealso \link[ggplot2]{ggplot}, \link{PS}

profileplot <- function(form,person.id,standardize=TRUE,interval=10,by.pattern=TRUE,original.names=TRUE) {
	
	if(standardize) {form <- scale(form)}
	
	if(by.pattern) {
	
	subscore.id <- subscore <- NULL
	n <- ncol(form)
	k <- nrow(form)
	
	if(original.names){labels <- colnames(form)}
	else {labels <- c(paste("v",1:n,sep=""))}
	
	#Level scores
	level <- as.matrix(rowMeans(form,na.rm = TRUE),ncol=1,nrow=k)
	
	#Pattern scores
	pattern <- matrix(,ncol=n, nrow=k)
	for (i in 1:n) {
		pattern[,i] <- form[,i]-level
	}
	
	if (missing(person.id)) {
		id=c(1:k)
		form1 <- cbind(id,form)
		level1 <- cbind(id,level)
		pattern1 <- cbind(id,pattern)
		
		colnames(form1) <- c("id",paste("s",c(1:n),sep=""))
		colnames(level1) <- c("id","level")
		colnames(pattern1) <- c("id",paste("p",c(1:n),sep=""))
		
		form1 <- as.data.frame(form1)
		pattern1 <- as.data.frame(pattern1)
		level1 <- as.data.frame(level1)
		
		form2 <- melt(form1, id="id")
		pattern2 <- melt(pattern1, id="id")
		colnames(form2) <- c("id","subscore.id","subscore")
		colnames(pattern2) <- c("id","pattern.id","pattern")
		
		form.long <- merge(form2, pattern2, by="id", all=TRUE)
		form.long <- merge(form.long, level1, by="id", all=TRUE)
		form.long$id <- as.factor(form.long$id)
		
	} 
	else {
		id=as.data.frame(person.id)
		form1 <- cbind(id,form)
		level1 <- cbind(id,level)
		pattern1 <- cbind(id,pattern)
		
		colnames(form1) <- c("id",paste("s",c(1:n),sep=""))
		colnames(level1) <- c("id","level")
		colnames(pattern1) <- c("id",paste("p",c(1:n),sep=""))
		
		form1 <- as.data.frame(form1)
		pattern1 <- as.data.frame(pattern1)
		level1 <- as.data.frame(level1)
		
		form2 <- melt(form1, id="id")
		pattern2 <- melt(pattern1, id="id")
		colnames(form2) <- c("id","subscore.id","subscore")
		colnames(pattern2) <- c("id","pattern.id","pattern")
		
		form.long <- merge(form2, pattern2, by="id", all=TRUE)
		form.long <- merge(form.long, level1, by="id", all=TRUE)
		form.long$id <- as.factor(form.long$id)
	}
	
	
	max.s <- max(form.long$subscore, na.rm = TRUE)
	min.s <- min(form.long$subscore, na.rm = TRUE)
	int <- (max.s-min.s)/interval
	
	plot1 <- ggplot(form.long, aes(x=subscore.id, y=subscore, group=id, color=id)) + geom_line() + geom_point(size=3, fill="white") + scale_colour_hue(name="Person",l=30) +
		theme(panel.background = element_rect(fill='white', colour='black'),panel.grid.minor=element_blank(), 
					panel.grid.major=element_blank()) + scale_x_discrete(name=" ", labels=labels)+
		scale_y_continuous(name="Scores", limits=c(min.s,max.s),breaks=seq(min.s, max.s, int)) + scale_shape_discrete(name="Person")}
	
	else {
	
  
  	
  	#Prepare data for plotting
  	form <- as.data.frame(form)
  	numvariables <- ncol(form)
  	colours <- brewer.pal(numvariables,"Set1")
  	mymin <- 1e+20
  	mymax <- 1e-20
  
  	for (i in 1:numvariables){
    		Scores <- form[,i]
		mini <- min(Scores)
    		maxi <- max(Scores)
    		if (mini < mymin) { mymin <- mini }
    		if (maxi > mymax) { mymax <- maxi }
  	}
  
  	if(original.names) {names <- colnames(form)}
  	else {names <- c(paste("v",1:numvariables,sep=""))}
  
  
  	#Plot the variables
  	for (i in 1:numvariables) {
		Scoresi <- form[,i]
    		namei <- names[i]
    		colouri <- colours[i]
    		if (i == 1) { plot1 <- plot(Scoresi,col=colouri,type="l",ylim=c(mymin,mymax),ylab="Score",xlab="Index") }
    		else         {points(Scoresi, col=colouri,type="l")                                     }
    		lastxval <- length(Scoresi)
    		lastyval <- Scoresi[length(Scoresi)]
    		text((lastxval),(lastyval),namei,col="black",cex=0.6)}
	}
    return(plot1)
}

