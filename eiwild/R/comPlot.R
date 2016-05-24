#' Plot Diagnostics to compare eiwild objects
#' 
#' @description
#' This function uses plot diagnostics (see \code{\link[eiwild]{plot.eiwild}}) to compare different
#' eiwild objects using functions of coda package
#' 
#' @param eiList list of \code{eiwild} objects
#' @param whichCell which cell to plot
#' @param whichPlot which type of plot. see Details for more information
#' @param whichParam which parameter should be plotted \code{"alphaDraws"} or \code{"cellCounts"}
#' @param rollLim specifying other \code{ylim}-values for rolling mean plot (Default=\code{NULL})
#' @param rollCol specifying other \code{col}-values for rolling mean plot (Default=\code{NULL})
#' @param ... arguments given to corresponding coda function
#' 
#' @details
#' \code{whichPlot} controls the plot diagnostic to run:
#' \itemize{
#'  \item \code{1} passes arguments to \code{\link[coda]{traceplot}}
#'  \item \code{2} passes arguments to \code{\link[lattice]{densityplot}}
#'  \item \code{3} calculates Running Mean with \code{eiwild:::rollMean}
#'  \item \code{4} passes arguments to \code{\link[coda]{gelman.plot}}.
#'           Output of \code{\link[coda]{gelman.diag}} will be title of this plot.
#' }
#' 
#' @seealso
#' \code{\link[coda]{mcmc}}
#' \code{\link[eiwild]{plot.eiwild}}
#' \code{\link[eiwild]{indAggEi}}
#' 
#' @examples
#' \dontrun{
#' # loading some fake election data
#' data(topleveldat)
#' form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
#' set.seed(1234)
#' out1 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                 sample=1000, thinning=2, burnin=100, verbose=100)
#'out2 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                 sample=1000, thinning=2, burnin=100, verbose=100)
#'out3 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                 sample=1000, thinning=2, burnin=100, verbose=100)
#'out4 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                 sample=1000, thinning=2, burnin=100, verbose=100)
#' 
#' eiList <- list(out1, out2, out3, out4)
#' 
#' comPlot(eiList, whichCell=1, whichPlot=1)
#' comPlot(eiList, whichCell="counts.CSU_1.CSU_2", whichPlot=1)
#' comPlot(eiList, whichCell=1, whichPlot=1, smooth=TRUE)
#' 
#' comPlot(eiList, whichCell=1, whichPlot=2)
#' 
#' comPlot(eiList, whichCell=1, whichPlot=3)
#' 
#' comPlot(eiList, whichCell=1, whichPlot=4)
#' comPlot(eiList, whichCell=1, whichPlot=4)
#' comPlot(eiList, 1, 3, whichParam="alphaDraws")
#' 
#' comPlot(eiList, "alpha.CSU_1.CSU_2", 3, whichParam="alphaDraws")
#' }
#' 
#' @export


comPlot <- function(eiList, whichCell, whichPlot, whichParam="cellCounts",
                    rollLim=NULL, rollCol=NULL, ...){
   
    ## checking stuff
    if( class(eiList)[1]!= "list")
        stop("\"eiList\" has to be list-object!", call.=FALSE)
    if(length(eiList)>1){
        le <- length(eiList)
    } else { 
        stop("\"eiList\" has to be of length bigger than 1. If not use standard plot-function", call. = FALSE) 
    }
    
    if(!whichParam %in% c("cellCounts", "alphaDraws"))
        stop("whichParam has to be \"alphaDraws\" or \"cellCounts\"!", call. = FALSE)
    
    if(all(sapply(1:length(eiList), function(j) class(eiList[[j]])[1]) %in% "eiwild")){
        objList <- vector("list",le)
        for(i in 1:le)
            objList[[i]] <- eiList[[i]]$draws[[whichParam]]
    } else { 
        stop("only \"eiwild\"-objects are allowed!", call. = FALSE)
    }
    
    for(i in 1:le){
        if(is.numeric(whichCell) && whichCell > ncol(objList[[i]])){
            stop("\"whichCell\" is not a Cell in eiwild-object!", call. = FALSE)
        } else if(is.character(whichCell) && !whichCell %in% colnames(objList[[i]])){
            stop("\"whichCell\" is not a Cell in eiwild-object!", call. = FALSE)
        }
    }
    
    if(is.numeric(whichCell)){ ## for title purposes
        varname <- colnames(objList[[1]])[whichCell] 
    } else{
        varname <- whichCell
    }
    
    ## making mcmc.list
    mcList <- mcmc.list()
    for(i in 1:le)
        mcList[[i]] <- objList[[i]][,whichCell]
    
    ## plotting objects
    if(whichPlot==1){ #traceplot
        traceplot(mcList, main=paste("Trace of",varname),...)
        
    } else if(whichPlot==2){#densityplot
        densityplot(mcList, main=paste("Density of",varname), ...)
        
    } else if(whichPlot==3){ #rolling mean
        means <- lapply(1:le, function(j) rollMean(mcList[[j]]))
        if(is.null(rollCol)){
            colors <- 1:le
        } else{
            colors <- rollCol
        }
        if(is.null(rollLim)){
            lim <- range(means)
        } else {
            lim <- rollLim
        }
        
        plot(1:length(means[[1]]), means[[1]],type="l",
             ylab="Means", xlab="Iterations", main=paste("Rolling mean of\n",varname),
             ylim=lim, col=colors[1], ...)          
        for(i in 2:le)
            lines(1:length(means[[i]]), means[[i]], col=colors[i])
    } else if(whichPlot==4){
        gelDiag <- gelman.diag(mcList)
        gelman.plot(mcList, main=paste(colnames(gelDiag$psrf),round(gelDiag$psrf,4),collapse="\n"),...)
    } else{
        stop("\"whichPlot\" has to be 1, 2, 3 or 4!", call. = FALSE)
    }
}


