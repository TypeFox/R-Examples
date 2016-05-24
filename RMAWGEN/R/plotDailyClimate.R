
NULL


#'   
#' Plots daily climatology through one year
#'  
#' @param data matrix whose columns contain daily-averaged climatic series of variables (e.g. maximum or minum daily averaged temperature obtained by spline interpolation of monthly climatology) 
#' @param origin origin date corresponding to the first row of \code{data}  
#' @param when   start day for daily climatology plot
#' @param title,xlab,ylab,col,lwd see \code{\link{plot.default}}
#' @param nday number of days in one year. Default is 365.
#' @param bicolor logical variable. If \code{TRUE} and \code{data} represents climatologies of minimun and maximum daily temperature, the lines are plotted with blue and red colors respectively.

#' 
#' @export 
#'   
#' @author  Emanuele Cordano, Emanuele Eccel
#'       
#' @return  a matrix containing the plotted variables


# THIS FUNCTION IS OBSOLATE!!!



plotDailyClimate <-
function (data,title="Daily_Avereged_Temperture_in_one_year",origin="1961-1-1",when="1979-1-1",ylab="Temperature [degC]",xlab="Time [days]",nday=365,bicolor=FALSE,col="black",lwd=1) { 
	

	
	
	
	
	data0  <- as.data.frame(extractdays(data,when=when,nday=nday,origin=origin))
	
	nc=ncol(data0)
	
	ylim <- c(min(data0),max(data0))
	
	if (bicolor) {
		colv=array("red",ncol(data0))
	} else {
		colv=array(col,ncol(data0))
	} 

	
	plot(1:nday,data0[,1],ylim=ylim,ylab=ylab,xlab=xlab,main=title,type="l",col=colv[1],lwd=lwd)
	
	
	for (i in 2:nc) {
		
		if ((bicolor) && (i>nc/2)) { 
			
			colv[i]="blue"
			
		} 
		
	
		
		lines(data0[,i],lty=i,col=colv[i],lwd=lwd)
		
		
	}
	
	legend("topleft",lty=1:nc,lwd=lwd,col=colv,legend=names(data0), inset = .05)
	
	
	return(data0)
	
	
}

