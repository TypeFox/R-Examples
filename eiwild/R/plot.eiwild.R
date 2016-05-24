#' Plot diagnostics for an eiwild Object
#' 
#' @description
#' plot diagnostics for an eiwild Object containing autocorrelation plots, density plots and trace
#' using functions of coda package
#' 
#' @param x an object of class \code{eiwild}
#' @param whichPlot which type of plot. see Details for more information
#' @param whichParam which parameter should be plotted \code{"alphaDraws"} or \code{"cellCounts"}
#' @param whichCell which cell to plot
#' @param layout logical if automatic layout of plot should be made with help of
#'           \code{par=mfrow()}
#' @param ... further graphical parameters given to correspondent coda function
#' 
#' @details 
#' \code{plot.eiwild} uses the plot diagnostic functions of the coda package. 
#' The default is \code{NULL} which passes arguments to \code{plot.mcmc}:
#' \itemize{
#'  \item \code{1} passes arguments to \code{\link[coda]{traceplot}}
#'  \item \code{2} passes arguments to \code{\link[coda]{autocorr.plot}}
#'  \item \code{3} passes arguments to \code{\link[coda]{densplot}}
#'  \item \code{4} calculates and plots rolling mean
#' }
#' 
#' @seealso
#' \code{\link[coda]{plot.mcmc}}
#' \code{\link[coda]{mcmc}}
#' 
#' @examples
#' \dontrun{
#' # loading some fake election data
#' data(topleveldat)
#' form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
#' set.seed(1234)
#' res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                 sample=1000, thinning=2, burnin=100,verbose=100)
#' 
#' plot(res, whichPlot=1)
#' plot(res, whichPlot=2)
#' plot(res, whichPlot=3)
#' plot(res, whichPlot=4)
#' 
#' plot(res, whichPlot=1, whichCell=c(1,4,5))
#' plot(res, whichPlot=2, whichCell=c(1,4,5))
#' plot(res, whichPlot=3, whichCell=c(1,4,5))
#' plot(res, whichPlot=4, whichCell=c(1,4,5))
#' 
#' plot(res, whichPlot=1, whichCell=c(1))
#' plot(res, whichPlot=2, whichCell=c(1))
#' plot(res, whichPlot=3, whichCell=c(1))
#' plot(res, whichPlot=4, whichCell=c(1))
#' 
#' plot(res, whichPlot=1, whichParam="alphaDraws")
#' plot(res, whichPlot=2, whichParam="alphaDraws")
#' plot(res, whichPlot=3, whichParam="alphaDraws")
#' plot(res, whichPlot=4, whichParam="alphaDraws")
#' 
#' par(mfrow=c(2,2))
#' plot(res, whichPlot=1, whichCell=1, layout=FALSE)
#' plot(res, whichPlot=2, whichCell=1, layout=FALSE)
#' plot(res, whichPlot=3, whichCell=1, layout=FALSE)
#' plot(res, whichPlot=4, whichCell=1, layout=FALSE)
#' }
#' 
#' @method plot eiwild
#' @export

plot.eiwild <- function(x, whichPlot=NULL, whichParam="cellCounts", whichCell=NULL, 
                        layout=TRUE, ...){
    
    if(whichParam %in% c("alphaDraws","cellCounts")){
        obj <- x$draws[[whichParam]]
    } else { stop("whichParam has to be \"alphaDraws\" or \"cellCounts\"!") }
    
    if(is.null(whichCell)){
        ##all cells
        cells <- 1:ncol(obj) 
    } else if(is.numeric(whichCell)){
        if(all(whichCell %in% 1:ncol(obj))){
            ## specific cells
            cells <- whichCell 
        } else { 
            stop(paste("whichCell has to be numbers within range",1,"to",ncol(obj)))
        }
    } else { 
        stop("whichCell has to be NULL (default) or numeric vector") 
    }
    
    iter <- nrow(obj)
    if(length(cells)==1){ ## because titles and obj[,i] doesn't work for one parameter
        temp <- t(t(obj[,cells]))
        colnames(temp) <- colnames(obj)[cells]  
        obj <- temp
    } else {
        obj <- obj[,cells]
    }
    
    tit <- colnames(obj)
    
    if(is.null(whichPlot)){
        ### Default plot function
        plot(obj,...)
        
    } else if(whichPlot==1) {
        ### Traceplot
        if(layout==TRUE){
        par(mfrow=calcmfrow(length(cells)))
        for(i in 1:length(cells))
            traceplot(obj[,i], main=paste("Trace of",tit[i]),...)
        par(mfrow=c(1,1))
        } else {
            traceplot(obj, main=paste("Trace of",tit),...)
        }
        
    } else if(whichPlot==2){
        ### autcorrplot
        if(layout==TRUE){
            par(mfrow=calcmfrow(length(cells)))
            for(i in 1:length(cells))
                autocorr.plot(obj[,i], auto.layout=FALSE, main=tit[i],...)
            par(mfrow=c(1,1))
        } else {
            autocorr.plot(obj, auto.layout=FALSE, main=tit,...)
        }
        
    } else if(whichPlot==3){
        ### densplot
        if(layout==TRUE){
        par(mfrow=calcmfrow(length(cells)))
        for(i in 1:length(cells))
            densplot(obj[,i], main=paste("Density of",tit[i]),...)
        par(mfrow=c(1,1))
        } else {
            densplot(obj, main=paste("Density of",tit),...)
        }
        
    } else if(whichPlot==4){
        ### Rolling Mean
        tit <- paste("Rolling mean of\n",colnames(obj))
        if(layout==TRUE){
        par(mfrow=calcmfrow(length(cells)))
        for(i in 1:length(cells))
            plot(1:iter, rollMean(obj[,i]), type="l", main=tit[i],
                 xlab="Iterations", ylab="means",...)
        par(mfrow=c(1,1))
        } else {
            plot(1:iter, rollMean(obj), type="l", main=tit,
                 xlab="Iterations", ylab="means",...)
        }
        
    } else {
        stop("whichPlot has to be NULL, 1, 2, 3 or 4!")
    }
}


