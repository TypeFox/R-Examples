# file plot_RMAWGENts.R
# 
# This file contains function to generate Q-Q plots of   stochastically generated daily temperature vs daily observed temperature. 
#
#
# author: Emanuele Cordano on 12-01-2012
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################
NULL
#'
#' It makes the Q-Q plots observed vs generated time series of daily maximum, minimum temperature and daily thermal range for a list of collected stochastic generations 
#' 
#' @param Tx_mes data frame containing measured daily maximum temperature    
#' @param Tn_mes data frame containing measured daily minimum temperature 
#' @param prec_mes data frame containing measured daily precipitation (in millimeters) 
#' @param Tx_spline data frame containing spline-interpolated daily maximum temperature. Default is \code{NULL} and not considered for Q-Q plot.    
#' @param Tn_spline data frame containing spline-interpolated daily minimum temperature  Default is \code{NULL} and not considered for Q-Q plot. 
#' @param Tx_gen data frame containing generated daily maximum temperature 
#' @param Tn_gen data frame containing generated daily minimum temperature 
#' @param prec_gen data frame containing generated daily precipitation (in millimeters) 
#' @param when day indices on which the data frame are extracted for Q-Q plot. Default is \code{1:nrow(Tn_mes)} (in \code{qqplot_RMAWGEN_Tn}) or \code{1:nrow(Tx_mes)} (otherwise) 
#' @param xlab,ylab lables of \code{x} and \code{y} axes. See \code{\link{qqplot}}.
#' @param station identification name (ID) of the station used for the Q-Q plot
#' @param main main titles for each plot. Default is \code{names(Tn_gen)}  (in \code{qqplot_RMAWGEN_Tn}) or \code{names(Tx_gen)} (otherwise)
#' @param pdf name of pdf file if output is written in a pdf file 
#' @param xlim see \code{\link{qqplot}}. Default is \code{range(Tn_mes)} (in \code{qqplot_RMAWGEN_Tn}) or \code{range(Tx_mes)} (in \code{qqplot_RMAWGEN_Tx})  .or \code{range(Tx_mes-Tn_mes)} (in \code{qqplot_RMAWGEN_deltaT})
#' @param ylim,cex,cex.main,cex.lab,cex.axis see \code{\link{qqplot}} and  \code{\link{plot}}
#' @param lag lag (current index included) on whose value  the precipitation addition is made. See \code{\link{qqplot.lagged}}.
#' 
#' 
#' @note \code{Tx_gen},{Tn_gen} and \code{main} must have an even number of elements.
#' 
#' @rdname qqplot_RMAWGEN_Tx
#' @export
#' @author Emanuele Cordano
#' 
#' 
#' 
#'  

qqplot_RMAWGEN_Tx <- function (Tx_mes,Tx_gen,Tn_gen,Tn_mes,Tx_spline=NULL,Tn_spline=NULL,xlab="observed",ylab="simulated",when=1:nrow(Tx_mes),main=names(Tx_gen),station,pdf=NULL,xlim=range(Tx_mes),ylim=xlim,cex=0.4,cex.main=1.0,cex.lab=1.0,cex.axis=1.0){
	
	if (!is.null(pdf)) pdf(pdf)
	
	N <- length(main)
	Q <- as.integer(N/2) 
	par(mfrow=c(Q,Q))
	for(i in 1:N) {
		
		if (is.null(Tx_spline))  {
			qqplot(Tx_mes[when,station],Tx_gen[[i]][when,station],xlab=xlab,ylab=ylab,main=main[i],cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim,ylim=ylim)
		} else {
			qqplot(Tx_mes[when,station]-Tx_spline[when,station],Tx_gen[[i]][when,station]-Tx_spline[when,station],xlab=xlab,ylab=ylab,main=main[i],cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim,ylim=ylim)
		}
			abline(0,1)
		
	}
	if (!is.null(pdf)) dev.off()
	
	
	
}

NULL
#'
#' @rdname qqplot_RMAWGEN_Tx
#' @export

qqplot_RMAWGEN_Tn  <- function (Tx_mes,Tx_gen,Tn_gen,Tn_mes,Tx_spline=NULL,Tn_spline=NULL,xlab="observed",ylab="simulated",when=1:nrow(Tn_mes),main=names(Tn_gen),station,pdf=NULL,xlim=range(Tn_mes),ylim=xlim,cex=0.4,cex.main=1.0,cex.lab=1.0,cex.axis=1.0){
	
	if (!is.null(pdf)) pdf(pdf)
	
	N <- length(main)
	Q <- as.integer(N/2) 
	par(mfrow=c(Q,Q))
	for(i in 1:N) {
		
		if (is.null(Tn_spline))  {
			qqplot(Tn_mes[when,station],Tn_gen[[i]][when,station],xlab=xlab,ylab=ylab,main=main[i],cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim,ylim=ylim)
		} else {
			qqplot(Tn_mes[when,station]-Tn_spline[when,station],Tn_gen[[i]][when,station]-Tn_spline[when,station],xlab=xlab,ylab=ylab,main=main[i],cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim,ylim=ylim)
		}
	abline(0,1)
		
	}
	if (!is.null(pdf)) dev.off()
	
	
	
}

NULL
#'
#' @rdname qqplot_RMAWGEN_Tx
#' @export

qqplot_RMAWGEN_deltaT  <- function (Tx_mes,Tx_gen,Tn_gen,Tn_mes,xlab="observed",ylab="simulated",when=1:nrow(Tx_mes),main=names(Tx_gen),station,pdf=NULL,xlim=range(Tx_mes-Tn_mes),ylim=xlim,cex=0.4,cex.main=1.0,cex.lab=1.0,cex.axis=1.0){
	
	if (!is.null(pdf)) pdf(pdf)
	
	N <- length(main)
	Q <- as.integer(N/2) 
	par(mfrow=c(Q,Q))
	for(i in 1:N) {
		
		qqplot(Tx_mes[when,station]-Tn_mes[when,station],Tx_gen[[i]][when,station]-Tn_gen[[i]][when,station],xlab=xlab,ylab=ylab,main=main[i],cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim,ylim=ylim)
		abline(0,1)
		
	}
	if (!is.null(pdf)) dev.off()
	
	
	
}

NULL
#'
#' @rdname qqplot_RMAWGEN_Tx
#' @export

qqplot_RMAWGEN_prec  <- function (prec_mes,prec_gen,xlab="observed",ylab="simulated",when=1:nrow(prec_mes),main=names(prec_gen),station,pdf=NULL,xlim=range(prec_mes),ylim=xlim,cex=0.4,cex.main=1.0,cex.lab=1.0,cex.axis=1.0,lag=1){
	
	if (!is.null(pdf)) pdf(pdf)
	
	N <- length(main)
	Q <- as.integer(N/2) 
	par(mfrow=c(Q,Q))
	for(i in 1:N) {
		
		qqplot.lagged(x=prec_mes[when,station],y=prec_gen[[i]][when,station],lag=lag,xlab=xlab,ylab=ylab,main=main[i],cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,xlim=xlim,ylim=ylim)
		abline(0,1)
		
	}
	if (!is.null(pdf)) dev.off()
	
	
	
}


