marimapred <-
function(TimeSeries,TimeSeriesCont=NULL,n.ahead=NULL,na.action=na.omit,xreg=NULL,newxreg=NULL,se.fit=FALSE,plot=FALSE,range.p=0.2,ylab=NULL,xlab=NULL,main=NULL){
    if(!is.null(TimeSeriesCont)) {
		if(!is.null(n.ahead)) {
			if(!is.null(xreg) & !is.null(newxreg)) {
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, MoreArgs = list(se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, ylab=colnames(TimeSeries), n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, MoreArgs = list(se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
			else{
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, n.ahead=n.ahead, MoreArgs = list(xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, ylab=colnames(TimeSeries), n.ahead=n.ahead, MoreArgs = list(xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
		}
		else{
			if(!is.null(xreg) & !is.null(newxreg)) {
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, xreg=xreg, newxreg=newxreg, MoreArgs = list(n.ahead=n.ahead, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, ylab=colnames(TimeSeries), xreg=xreg, newxreg=newxreg, MoreArgs = list(n.ahead=n.ahead, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
			else{
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, MoreArgs = list(n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, TimeSeriesCont, ylab=colnames(TimeSeries), MoreArgs = list(n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
		}
	}
    else{
        if(!is.null(n.ahead)) {
			if(!is.null(xreg) & !is.null(newxreg)) {
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, MoreArgs = list(timeseries.cont=TimeSeriesCont, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, ylab=colnames(TimeSeries), n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, MoreArgs = list(timeseries.cont=TimeSeriesCont, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
			else{
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, n.ahead=n.ahead, MoreArgs = list(timeseries.cont=TimeSeriesCont, xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, ylab=colnames(TimeSeries), n.ahead=n.ahead, MoreArgs = list(timeseries.cont=TimeSeriesCont, xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
		}
		else{
			if(!is.null(xreg) & !is.null(newxreg)) {
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, xreg=xreg, newxreg=newxreg, MoreArgs = list(timeseries.cont=TimeSeriesCont, n.ahead=n.ahead, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, ylab=colnames(TimeSeries), xreg=xreg, newxreg=newxreg, MoreArgs = list(timeseries.cont=TimeSeriesCont, n.ahead=n.ahead, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
			else{
				if(!is.null(ylab)) {
					Predictions <- mapply(arimapred, TimeSeries, MoreArgs = list(timeseries.cont=TimeSeriesCont, n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, ylab=ylab, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
				else {
					Predictions <- mapply(arimapred, TimeSeries, ylab=colnames(TimeSeries), MoreArgs = list(timeseries.cont=TimeSeriesCont, n.ahead=n.ahead, xreg=xreg, newxreg=newxreg, se.fit=se.fit, na.action=na.action, range.p=range.p, plot=plot, xlab=xlab, main=main), SIMPLIFY = TRUE, USE.NAMES = TRUE)
				}
			}
		}
    }
    return (Predictions)
}