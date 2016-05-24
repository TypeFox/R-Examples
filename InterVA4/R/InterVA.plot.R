#' Summarize population level cause-specific mortality fraction as InterVA4 suggested.
#' 
#' The function takes input of a list of va object and calculates the 
#' cause-specific mortality fraction. It only calculates CSMF as aggregation of up to the third largest causes.
#' 
#' @param va The list of va object to summarize.
#' @return \item{dist.cod}{The cause-specific mortality fraction (including undetermined category).}
#' @author Zehang LI, Tyler McCormick, Sam Clark
#' @keywords interVA
#' @seealso \code{\link{CSMF}}
#' @examples
#' 
#' data(SampleInput)
#' sample.output <- InterVA(SampleInput, HIV = "h", Malaria = "v", directory = "VA test",
#'        filename = "VA_result", output = "extended", append = FALSE)
#' ## Get CSMF without plots
#' csmf<- CSMF.interVA4(sample.output$VA)
#' data(SampleInput)
#' 

CSMF.interVA4 <- function(va){
   # for future compatibility with non-standard input
    for(i in 1:length(va)){
        if(!is.null(va[[i]]$wholeprob)){
            causenames <- names(va[[i]]$wholeprob)
            causeindex <- 1:length(causenames)
            break
        }
    }
   
    # fix for removing the first 3 preg related death in standard input
    if(causenames[1] == "Not pregnant or recently delivered" &&
        causenames[2] == "Pregnancy ended within 6 weeks of death" &&
        causenames[3] == "Pregnant at death"){
            causeindex <- causeindex[-c(1:3)]
            causenames <- causenames[-c(1:3)]    
    }


    ## Check if there is a valid va object
    if(length(va) < 1){
        cat("No va object found")
        return
    }
    ## Initialize the population distribution
    dist <- NULL
    for(i in 1:length(va)){
        if(!is.null(va[[i]][14])){  
            dist <- rep(0, length(unlist(va[[i]][14]))) 
            break
        }
    } 
    undeter <- 0

    ## pick not simply the top 3 causes, but the top 3 causes reported by InterVA
    for(i in 1:length(va)){
        if(is.null(va[[i]][14])) next
        this.dist <- unlist(va[[i]][14])
        if(max(this.dist) < 0.4){
          undeter <- undeter + sum(this.dist)  
        }else{
            cutoff.3 <- this.dist[order(this.dist, decreasing = TRUE)[3]]
            cutoff.2 <- this.dist[order(this.dist, decreasing = TRUE)[2]]
            cutoff.1 <- this.dist[order(this.dist, decreasing = TRUE)[1]]
            cutoff <- min(max(cutoff.1 * 0.5 , cutoff.3), max(cutoff.1 * 0.5 , cutoff.2))

            undeter <- undeter + sum(this.dist[which(this.dist < cutoff)])
            this.dist[which(this.dist < cutoff)] <- 0
            if(!is.null(va[[i]][14])) dist <- dist + this.dist
        }
    } 
    ## Normalize the probability for CODs
    if(undeter > 0){
        dist.cod <- c(dist[causeindex], undeter)
        dist.cod <- dist.cod/sum(dist.cod)
        names(dist.cod)<-c(causenames, "Undetermined")
    }else{
        dist.cod <- dist[causeindex]/sum(dist[causeindex])
        names(dist.cod)<-causenames
    } 

    return(dist.cod)  
}

#' Summarize and plot a population level distribution of va probabilities.
#' 
#' This function has been deprecated as of version 1.6. Use 'CSMF' instead.
#' 
#' 
#' @param va The list of va object to summarize.
#' @param top.aggregate Integer indicating how many causes from the top need to go into
#' summary. The rest of the probabilities goes into an extra category
#' "Undetermined".  When set to NULL, default is all causes to be considered.
#' This is only used when \code{InterVA.rule} set to "FALSE".
#' @param InterVA.rule If it is set to "TRUE", only the top 3 causes reported by 
#' InterVA4 is calculated into CSMF as in InterVA4. The rest of probabilities 
#' goes into an extra category "Undetermined". Default set to "FALSE".
#' @param noplot A logical value indicating whether the plot will be shown. If
#' it is set to "TRUE", only the CSMF will be returned.
#' @param min.prob The minimum probability that is to be plotted in bar chart,
#' or to be labeled in pie chart.
#' @param type An indicator of the type of chart to plot.  "pie" for pie chart;
#' "bar" for bar chart.
#' @param ... Arguments to be passed to/from graphic function
#' \code{\link[graphics]{barplot}}, \code{\link[graphics]{pie}}, and more
#' graphical paramters (see \code{\link[graphics]{par}}). They will affect the
#' main title, size and font of labels, and the radius of the pie chart.
#' @return \item{dist.cod}{The population probability of CODs.}
#' @author Zehang LI, Tyler McCormick, Sam Clark
#' @seealso \code{\link{CSMF.interVA4}}
#' @keywords interVA
#' 
#' 
Population.summary<-function (va, top.aggregate = NULL, InterVA.rule = FALSE, noplot = FALSE, type="bar",  min.prob = 0.01, ... ) {
    .Deprecated("CSMF")
    CSMF(va, top.aggregate = NULL, InterVA.rule = FALSE, noplot = FALSE, type="bar",  min.prob = 0.01, ... )
}


#' Summarize and plot a population level distribution of va probabilities.
#' 
#' The function takes input of a list of va object and produces a summary plot
#' for the population distribution.
#' 
#' 
#' @param va The list of va object to summarize.
#' @param top.aggregate Integer indicating how many causes from the top need to go into
#' summary. The rest of the probabilities goes into an extra category
#' "Undetermined".  When set to NULL, default is all causes to be considered.
#' This is only used when \code{InterVA.rule} set to "FALSE".
#' @param InterVA.rule If it is set to "TRUE", only the top 3 causes reported by 
#' InterVA4 is calculated into CSMF as in InterVA4. The rest of probabilities 
#' goes into an extra category "Undetermined". Default set to "FALSE".
#' @param noplot A logical value indicating whether the plot will be shown. If
#' it is set to "TRUE", only the CSMF will be returned.
#' @param top.plot the maximum number of causes to plot in bar plot
#' @param min.prob The minimum probability that is to be plotted in bar chart,
#' or to be labeled in pie chart.
#' @param type An indicator of the type of chart to plot.  "pie" for pie chart;
#' "bar" for bar chart.
#' @param ... Arguments to be passed to/from graphic function
#' \code{\link[graphics]{barplot}}, \code{\link[graphics]{pie}}, and more
#' graphical paramters (see \code{\link[graphics]{par}}). They will affect the
#' main title, size and font of labels, and the radius of the pie chart.
#' @return \item{dist.cod}{The population probability of CODs.}
#' @author Zehang LI, Tyler McCormick, Sam Clark
#' @seealso \code{\link{CSMF.interVA4}}
#' @keywords interVA
#' @examples
#' 
#' data(SampleInput)
#' sample.output <- InterVA(SampleInput, HIV = "h", Malaria = "v", directory = "VA test", 
#'                          filename = "VA_result", output = "extended", append = FALSE)
#' 
#' ## Get CSMF without plots
#' population.summary <- CSMF(sample.output$VA, noplot = TRUE)
#' 
#' 
#' ## Get CSMF by considering only top 3 causes for each death.
#' population.summary <- CSMF(sample.output$VA, top.aggregate = 3, noplot = TRUE)
#' 
#' ## Get CSMF by considering only top 3 causes reported by InterVA.  
#' ## This is equivalent to using CSMF.interVA4() command Note that
#' ## it's different from using all top 3 causses, since they may not 
#' ## all be reported 
#' CSMF.summary <- CSMF(sample.output, InterVA.rule = TRUE, 
#'    noplot = TRUE)
#' 
#' ## Population level summary using pie chart
#' CSMF.summary2 <- CSMF(sample.output, type = "pie", 
#'  min.prob = 0.01, main = "population COD distribution using pie chart", 
#'  clockwise = FALSE, radius = 0.7, cex = 0.7, cex.main = 0.8)
#' 
#' ## Population level summary using bar chart
#' CSMF.summary3 <- CSMF(sample.output, type = "bar", 
#'   min.prob = 0.01, main = "population COD distribution using bar chart", 
#'   cex.main = 1)
## Population level summary specified by number of top causes
#' CSMF.summary4 <- CSMF(sample.output, type = "bar", 
#'   top.plot = 5, main = "Top 5 population COD distribution", 
#'   cex.main = 1)
#' 
CSMF <-function (va, top.aggregate = NULL, InterVA.rule = FALSE, noplot = FALSE, type="bar",  top.plot = 10, min.prob = 0, ... ) {
	
    ## Check if there is a valid va object
    if(class(va) == "interVA"){
        va <- va$VA
    }

    # for future compatibility with non-standard input
    for(i in 1:length(va)){
        if(!is.null(va[[i]]$wholeprob)){
            causenames <- names(va[[i]]$wholeprob)
            causeindex <- 1:length(causenames)
            break
        }
    }
    
    # fix for removing the first 3 preg related death in standard input
    if(causenames[1] == "Not pregnant or recently delivered" &&
        causenames[2] == "Pregnancy ended within 6 weeks of death" &&
        causenames[3] == "Pregnant at death"){
            causeindex <- causeindex[-c(1:3)]
            causenames <- causenames[-c(1:3)]    
    }

    if(length(va) < 1){
		cat("No va object found")
		return
	}
    ## Initialize the population distribution
    dist <- NULL
    for(i in 1:length(va)){
        if(!is.null(va[[i]][14])){	
	        dist <- rep(0, length(unlist(va[[i]][14])))	
	        break
        }
    } 
    ## determine how many causes from top need to be summarized
    if(is.null(top.aggregate)) top.aggregate <- length(causeindex)
    undeter <- 0

    if(is.null(dist)){cat("No va probability found in input"); return}   
    ## Add the probabilities together
	if(!InterVA.rule){
        for(i in 1:length(va)){
            if(is.null(va[[i]][14])) {undeter = undeter + 1; next}
            this.dist <- unlist(va[[i]][14])
            cutoff <- this.dist[order(this.dist, decreasing = TRUE)[top.aggregate]]
            undeter <- undeter + sum(this.dist[which(this.dist < cutoff)])
            
            this.dist[which(this.dist < cutoff)] <- 0
            if(!is.null(va[[i]][14])) dist <- dist + this.dist
        }  
            ## Normalize the probability for CODs
        if(undeter > 0){
            dist.cod <- c(dist[causeindex], undeter)
            dist.cod <- dist.cod/sum(dist.cod)
            names(dist.cod)<-c(causenames, "Undetermined")
        }else{
            dist.cod <- dist[causeindex]/sum(dist[causeindex])
            names(dist.cod)<-causenames
        }      
    }else{
        dist.cod <- CSMF.interVA4(va)   
    }


    ## Check if there is CODs above the minimum cut-off for prob
    if(max(dist.cod) < min.prob){
        cat("No COD larger than the minimum probability cut off line")
        return
    }
    if(noplot){
    		return(dist.cod)
    }

    if(!is.null(top.plot)){
        if(top.plot < length(dist.cod)){
            thre <- sort(dist.cod, decreasing=TRUE)[top.plot]
            min.prob <- max(min.prob, thre)
        }
    }

    ## Make pie plot upon request
    if( type == "pie" ){
        dist.cod.sort <- sort(dist.cod, decreasing=TRUE)
        pie.color <- grey.colors(length(dist.cod.sort[dist.cod.sort >= min.prob]))
        pie.color.left <- rep(pie.color[length(pie.color)], length(dist.cod.sort[dist.cod.sort < min.prob]))
        pie.color <- c(pie.color, pie.color.left)
        pie(dist.cod.sort, col = pie.color,labels = names(dist.cod.sort)[dist.cod.sort > min.prob], ...)
        
    }
    ## Make bar plot upon request
    if( type == "bar"){
        dist.cod.min <- dist.cod[dist.cod >= min.prob ]
        dist.cod.min <- sort(dist.cod.min, decreasing = FALSE)
        par(las = 2)
        par(mar = c(5,15,4,2))
        bar.color <- grey.colors(length(dist.cod.min))
        bar.color <- rev(bar.color)
        barplot(dist.cod.min , horiz = TRUE,names.arg = names(dist.cod.min), col = bar.color, cex.names=0.8, xlab = "Probability", ...)
    }
    ## Save the population distribution
    dist.cod
	
}






#' Plot a individual level distribution of va probabilities.
#' 
#' The function takes input of a single va object and produces a summary plot
#' for it.
#' 
#' 
#' @param va A va object
#' @param min.prob The minimum probability that is to be plotted in bar chart,
#' or to be labeled in pie chart.
#' @param type An indicator of the type of chart to plot.  "pie" for pie chart;
#' "bar" for bar chart.
#' @param ... Arguments to be passed to/from graphic function
#' \code{\link[graphics]{barplot}}, \code{\link[graphics]{pie}}, and more
#' graphical paramters (see \code{\link[graphics]{par}}). They will affect the
#' main title, size and font of labels, and the radius of the pie chart.
#' @seealso \code{\link{CSMF}}
#' @keywords InterVA
#' @examples
#' 
#' data(SampleInput)
#' sample.output <- InterVA(SampleInput, HIV = "h", Malaria = "v", directory = "VA test", 
#'     filename = "VA_result", output = "extended", append = FALSE)
#' 
#' ## Individual level summary using pie chart
#' InterVA.plot(sample.output$VA[[7]], type = "pie", min.prob = 0.01, 
#'     main = "1st sample VA analysis using pie chart", clockwise = FALSE, 
#'     radius = 0.6, cex = 0.6, cex.main = 0.8)
#' 
#' 
#' ## Individual level summary using bar chart
#' InterVA.plot(sample.output$VA[[7]], type = "bar", min.prob = 0.01, 
#'     main = "2nd sample VA analysis using bar chart", cex.main = 0.8)
#' 
InterVA.plot <- function(va, type="bar", min.prob = 0.01, ... ){
    
    # for future compatibility with non-standard input
    if(!is.null(va$wholeprob)){
        causenames <- names(va$wholeprob)
        causeindex <- 1:length(causenames)
    }else{
        cat("Cause of death undetermined for this case\n")
        return
    }

    # fix for removing the first 3 preg related death in standard input
    if(causenames[1] == "Not pregnant or recently delivered" &&
        causenames[2] == "Pregnancy ended within 6 weeks of death" &&
        causenames[3] == "Pregnant at death"){
            causeindex <- causeindex[-c(1:3)]
            causenames <- causenames[-c(1:3)]    
    }


    ## Check if there is a valid va object
	if(length(va) < 1){
		cat("No va object found")
		return
	}
    ## Find the probability distribution
	dist <- unlist(va[14])
    dist.cod <- dist[causeindex]/sum(dist[causeindex])
    ## Check if there is CODs above the minimum cut-off for prob
    if(max(dist.cod) < min.prob){
        cat("No COD larger than the minimum probability cut off line")
        return
    }
    names(dist.cod)<-causenames
    ## Make pie plot upon request    
    if( type == "pie" ){
        dist.cod.sort <- sort(dist.cod, decreasing=TRUE)
        pie.color <- grey.colors(length(dist.cod.sort[dist.cod.sort >= min.prob]))
        pie.color.left <- rep(pie.color[length(pie.color)], length(dist.cod.sort[dist.cod.sort < min.prob]))
        pie.color <- c(pie.color, pie.color.left)
        pie(dist.cod.sort, col = pie.color,labels = names(dist.cod.sort)[dist.cod.sort > min.prob], ...)
        
    }
    ## Make bar plot upon request
    if( type == "bar"){
        dist.cod.min <- dist.cod[dist.cod >= min.prob ]
        dist.cod.min <- sort(dist.cod.min, decreasing = FALSE)
        par(las = 2)
        par(mar = c(5,15,4,2))
        bar.color <- grey.colors(length(dist.cod.min))
        bar.color <- rev(bar.color)
        barplot(dist.cod.min , horiz = TRUE,names.arg = names(dist.cod.min), col = bar.color, cex.names=0.8, xlab = "Probability", ...)
    }

}
