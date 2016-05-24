#' Visualize the Results of an APSIM Simulation
#' 
#' Plot routine for \code{"apsim"} class objects
#' 
#' Similar to the \code{\link[stats:plot.lm]{plot}} for \code{\link{lm}} objects, \code{plot.apsim} will 
#' plot each response in the results of an APSIM simulation on its own \code{ggplot2} object.  If the
#' \code{one_plot} argument is set to \code{TRUE} then \code{\link[ggplot2:facet_wrap]{facet_wrap}} is used to plot
#' all of the responses on one screen.  Alternatively, one response can be plotted by setting \code{y} to 
#' one variable corresponding to a column in \code{x}.
#' 
#' @name plot.apsim
#' @param x data frame of class \code{"apsim"} including the results of an APSIM simulation
#' @param y variable to plot on y-axis; if left empty all variables will be plotted on separate plots
#' @param ask logical; if \code{TRUE}, the user is asked before each plot, see \code{\link{par}(ask=.)} 
#' @param one_plot logical; if \code{TRUE} all variables are plotted on one faceted plot
#' @param ... additional arguments passed to \code{\link[ggplot2:qplot]{qplot}}
#' @S3method plot apsim
#' @method plot apsim
#' @export
#' @examples
#' \dontrun{
#' apsimExe <-"C:/Program Files (x86)/Apsim75-r3008/Model/Apsim.exe"
#' apsimWd <- "~/APSIM"
#' toRun <- c("Centro.apsim", "Continuous Wheat.apsim")
#' results <- apsim(exe = apsimExe, wd = apsimWd, files = toRun)
#' 
#' #Look at all of the results as a function of time in seperate plots
#' plot(results[[2]])
#' 
#' #Put all variables on one faceted plot
#' plot(results[[2]], one_plot = TRUE) + theme_bw()
#' 
#' #Plot just yield as a function of time
#' plot(results[[2]], y = 'yield') + geom_line(colour = 'red') + theme_bw()
#' }

plot.apsim<-function(x, y = NULL, ask = TRUE, one_plot = FALSE, ...){
  ncol<-ncol(x)
  cNames<-colnames(x)
  idNum<-which('Date'%in%cNames)
  
  if(is.null(y)){
    y<-cNames[-idNum]
  }else{
    
    possibleY <- try(match.arg(arg=y,choices=cNames[-idNum]),silent=T)
    
    if(class(possibleY)=="try-error"){
      stop(paste(y,"is not an available response; choose from:",paste(cNames,collapse=', ')))
    }else{
      y <- possibleY
    }
  }
  
  if(length(y)==1){

    return(qplot(x[,'Date'],x[,y],xlab='Date',ylab=y,...))
    
  }else if(one_plot==TRUE){
    
    mX <- melt(x,id.vars='Date')
    
    return(qplot(mX$Date,mX$value,data=mX,ylab="Value",xlab='Date',...)+facet_wrap(~variable,scales='free'))
    
    
  }else{
    
    if(ask){
      oask <- devAskNewPage(TRUE)
    }
    
    for(i in y){
      print(qplot(x[,'Date'],x[,i],xlab='Date',ylab=i,...)+theme_bw())
    }
    oask <- devAskNewPage(FALSE)
  }
}