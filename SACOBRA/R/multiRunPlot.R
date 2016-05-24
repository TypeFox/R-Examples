######################################################################################
# multiRunPlot
#
#' Plot the results from multipe COBRA runs.
#' 
#' Plot for each run one black curve 'error vs. iterations' and aggregate the mean curve (red) and 
#' the median curve (green) of all runs. 'error' is the distance between the ever-best feasible 
#' value and \code{optim}.
#' 
#' Print some diagnostic information: final median & mean error, percentage of runs which meet the
#' target (only if \code{optim} is available)). \cr 
#' 
#'  @param dfAll      the data frame of all runs, obtained with \code{\link{multiCOBRA}} or loaded 
#'                    from .Rdata file
#'  @param optim      [NULL] the true optimum (or best known value) of the problem (only for diagnostics).
#'                    If \code{optim==NULL}, we plot instead of errors the ever-best feasible values.
#'  @param fName      [""] the name of the .Rdata file, printed as subtitle
#'  @param main       [""] the name of the problem (e.g. "G01 problem"), printed as title  
#'  @param xlim       the x limits
#'  @param ylim       the y limits
#'  @param ylog       [TRUE]  logarithmic y-axis
#'  @param xlog       [FALSE] logarithmic x-axis
#'  @param target     [0.05] a single run meets the target, if the final error is smaller than \code{target}
#'  @param plotPDF    [FALSE] if TRUE, plot to 'fName'.pdf
#'  @param subPDF     [NULL] optional subdirectory where .pdf should go
#'  @param legendWhere  ["topright"]
#'
#'  @return \code{z3}, a vector containing for each run the ever-best feasible objective value
#'
#' @seealso   \code{\link{multiCOBRA}}, \code{\link{cobraPhaseII}}
#' @author Wolfgang Konen, Samineh Bagheri, Cologne Univeristy of Applied Sciences
#' @export
#' 
multiRunPlot <- function(dfAll,optim=NULL,fName="",main="",
                           xlim=NULL,ylim=c(1e-05,1e+04),ylog=TRUE, 
                           xlog=FALSE,target=0.05, plotPDF=FALSE, subPDF=NULL,
                           legendWhere="topright" ) 
{
  ffcMin = min(dfAll$ffc)
  ffcMax = max(dfAll$ffc)
  runList = unique(dfAll$run)
  everMax = max(dfAll$everBestFeas)
  if (is.null(optim)) {
    optimumAvail=FALSE
    optim=0.0
    ylab="best feasible"
  } else {
    optimumAvail=TRUE
    ylab="error"
  }
  if (is.null(xlim)) xlim = c(ffcMin,ffcMax)
  df = dfAll[dfAll$run==runList[1],]
  slog = ""
  if (xlog) slog=paste(slog,"x",sep="")
  if (ylog) slog=paste(slog,"y",sep="")
  pdfname = paste("results.d",sub(".Rdata",".pdf",fName),sep="/")
  if (!is.null(subPDF)) pdfname=paste(subPDF,pdfname,sep="/")
  if (plotPDF) {
    grDevices::pdf(pdfname,width=7,height=5.24)
    graphics::par(cex.axis=1.3)
    graphics::par(cex.lab=1.3)
    graphics::par(cex.sub=0.9)
  } else {
    graphics::par(cex.axis=1.0)
    graphics::par(cex.lab=1.0)
    graphics::par(cex.sub=0.7)    
  }
  
  #browser()
  graphics::plot(df$ffc,df$everBestFeas-optim,log=slog,type="l",main=main,sub=fName
       ,xlim=xlim,ylim=ylim,ylab=ylab,xlab="iterations") 
  
  for (i in runList) {
    df = dfAll[dfAll$run==runList[i],]
 
    # set dfAll$everBestFeas of each run to NA until the first feasible iterate arrives 
    ind <- which(df$feas==TRUE)
    if (length(ind)>0) {
      ind=min(ind)
      if (ind>1 & ind!=Inf)
        df$everBestFeas[1:(ind-1)] <- NA;   
      dfAll$everBestFeas[dfAll$run==runList[i]] <- df$everBestFeas;
    } else {
      cat(sprintf("NOTE: No feasible solution in run %d\n",i))
      dfAll$everBestFeas[dfAll$run==runList[i]] <- NA;      
    }   

    graphics::lines(df$ffc,df$everBestFeas-optim)
  }
  
  graphics::legend(legendWhere,c("mean","median"),lwd=2,col=c("red","green"))
  #browser()
  indInfeas = which(dfAll$feas==FALSE)
  z = stats::aggregate(dfAll$everBestFeas,list(dfAll$run),min, na.rm=TRUE, na.action=NULL)
  indInf = which(is.infinite(z$x))
  if (length(indInf)>0) {
    if (length(indInf)==1)
      cat(sprintf("*** There is %3d infinity solution (never feasible) ***\n",length(indInf)))
    if (length(indInf)>1)
      cat(sprintf("*** There are %3d infinity solutions (never feasible) ***\n",length(indInf)))
    cat("For the non-Inf solutions:\n")
    cat(sprintf("Avg solution: %10.5f +- %7.5f,  true minimum: %10.5f\n", 
                mean(z$x[-indInf]),sd(z$x[-indInf]),optim))  #,"\n"
  } else {
    cat(sprintf("Avg solution: %10.5f +- %7.5f,  true minimum: %10.5f\n", 
                mean(z$x),sd(z$x),optim))  #,"\n"
    
  }
  
  
  z = stats::aggregate(dfAll$everBestFeas,list(dfAll$ffc),mean,na.rm=T)
  names(z) <- c("ffc","everBestFeas")
  z = z[which(z$ffc<=ffcMax),]            # strip off ffc==41, if there
  graphics::lines(z$ffc,z$everBestFeas-optim,col="red",lwd=2)
  cat("Final avg error (mean):  ", z$everBestFeas[nrow(z)]-optim,"\n")  
  
  z2 = stats::aggregate(dfAll$everBestFeas,list(dfAll$ffc),stats::median,na.rm=T)
  names(z2) <- c("ffc","everBestFeas")
  z2 = z2[which(z2$ffc<=ffcMax),]            # strip off ffc==41, if there
  graphics::lines(z2$ffc,z2$everBestFeas-optim,col="green",lwd=2)
  cat("Final avg error (median):", z2$everBestFeas[nrow(z)]-optim,"\n")  
  
  z3 = dfAll$everBestFeas[dfAll$ffc==ffcMax]
  if (optimumAvail) {
    pTarget = length(which(z3-optim < target))/length(z3)
    cat("Target (",target,") reached with probability: ",pTarget,"\n")
  }
  
  if (plotPDF) {
    grDevices::dev.off()
    cat("Plot saved to",pdfname,"\n")
  }
  return(z3)
}
