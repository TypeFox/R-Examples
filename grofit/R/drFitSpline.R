drFitSpline <-
function (conc, test, drID="undefined", control=grofit.control())
{
# /// check input parameters
if (is(control)!="grofit.control")
     stop("control must be of class grofit.control!")

# /// conversion to handle even data.frame inputs
test <- as.vector(as.numeric(as.matrix(test)))
conc <- as.vector(as.numeric(as.matrix(conc)))
 
if(is.vector(conc)==FALSE ||is.vector(test)==FALSE)
     stop("drFitSpline: dose or response data must be a vector !")

if (control$neg.nan.act==FALSE){
    # /// removing missing data
    missings <- is.na(conc)|is.na(test)|!is.numeric(conc)|!is.numeric(test)
    conc     <- conc[!missings]
    test     <- test[!missings]

    # /// remove negative values
    negs <- (conc<0)|(test<0)
    conc <- conc[!negs]
    test <- test[!negs]
}
else{
    if (sum(is.na(conc)|is.na(test))) stop("drFitSpline: NA values encountered. Program terminated")
    if ((sum((conc<0))>0)|(sum((test<0))>0)) stop("drFitSpline: Negative values encountered. Program terminated")
    if ((FALSE%in%is.numeric(conc))||(FALSE%in%is.numeric(test))) stop("drFitSpline: Non numeric values encountered. Program terminated")
}

# /// check length of data after removal of values
if (length(test)<6){
	warning("drFitSpline: There is not enough valid data. Must have at least 6 unique values!")
	drFitSpline <- list(raw.conc = conc, raw.test=test, drID=drID, fit.conc=NA, fit.test=NA, spline = NA, parameters=list("EC50"=NA, "yEC50"=NA, "EC50.orig"=NA, "yEC50.orig"=NA), fitFlag = FALSE, reliable=NULL, control=control)
	class(drFitSpline)<-"drFitSpline"
	return(drFitSpline)
}

if (length(test)<control$have.atleast){
	warning("drFitSpline: number of valid data points is below the number specified in 'have.atleast'. See grofit.control().")
	drFitSpline <- list(raw.conc = conc, raw.test=test, drID=drID, fit.conc=NA, fit.test=NA, spline = NA, parameters=list("EC50"=NA, "yEC50"=NA, "EC50.orig"=NA, "yEC50.orig"=NA), fitFlag = FALSE, reliable=NULL, control=control)
	class(drFitSpline)<-"drFitSpline"
	return(drFitSpline)
}

# /// transformation of data...
if (control$log.x.dr==TRUE)  conc.log <- log(conc+1)
if (control$log.y.dr==TRUE)  test.log <- log(test+1)
if (control$log.x.dr==TRUE) {conc.fit <- log(conc+1)} else {conc.fit <- conc}
if (control$log.y.dr==TRUE) {test.fit <- log(test+1)} else {test.fit <- test}

# /// spline fit
spltest <- NULL
fitFlag <- TRUE
try(spltest<-smooth.spline(conc.fit,test.fit,spar=control$smooth.dr))
if (is.null(spltest)==TRUE){
    cat("Spline could not be fitted to data!\n")
    fitFlag <- FALSE
    if (is.null(control$smooth.dr)==TRUE){
        cat("This might be caused by usage of smoothing parameter NULL\n")
    }
    stop("Error in drFitSpline")	
}

# /// estimating EC 50 values
conc.min  <- min(conc.fit)
conc.max  <- max(conc.fit)
c.pred    <- seq(conc.min,conc.max,length.out=1000)
ytest     <- predict(spltest,c.pred)
yEC.test  <- (max(ytest$y)-min(ytest$y))/2.0+min(ytest$y)
last.test <- max(ytest$y)

kec.test<-1
for(k in 1:(length(c.pred)-1)){
  d1 <- (ytest$y[k]-yEC.test)
  d2 <- (ytest$y[k+1]-yEC.test)
  if (((d1<=0)&&(d2>=0))|((d1>=0)&&(d2<=0))){
      kec.test<-k
      break
  }
}
EC.test<-c.pred[kec.test]

# /// print EC50 to screen
if (control$suppress.messages == FALSE){
	cat("\n\n=== Dose response curve estimation ================\n")
	cat("--- EC 50 -----------------------------------------\n")
        cat(paste("-->", as.character(drID)))
	cat("\n")
	cat(paste(c("xEC50", "yEC50"),c(EC.test, yEC.test)))
}
	
if ((control$log.x.dr==TRUE)&&(control$log.y.dr==FALSE)){
	if (control$suppress.messages == FALSE){
		cat("\n--> Original scale \n")
		cat(paste(c("xEC50", "yEC50"),c(exp(EC.test)-1,yEC.test)))
	}
	EC.orig <- c(exp(EC.test)-1,yEC.test)
}
else{
	if ((control$log.x.dr==FALSE)&&(control$log.y.dr==TRUE)){
		if (control$suppress.messages == FALSE){
			cat("\n--> Original scale \n")
			cat(paste(c("xEC50", "yEC50"),c(EC.test,exp(yEC.test)-1)))
		}
	EC.orig <- c(EC.test,exp(yEC.test)-1)
	}
	else{
		if ((control$log.x.dr==TRUE)&&(control$log.y.dr==TRUE)){
			if (control$suppress.messages == FALSE){
				cat("\n--> Original scale \n")
				cat(paste(c("xEC50", "yEC50"),c(exp(EC.test)-1, exp(yEC.test)-1)))
			}
		EC.orig <- c(exp(EC.test)-1, exp(yEC.test)-1)
		}
		else{
			if ((control$log.x.dr==FALSE)&&(control$log.y.dr==FALSE)){
				EC.orig <- c(EC.test, yEC.test)
			}
		}
	}
}
if (control$suppress.messages == FALSE){ cat("\n\n\n") }

drFitSpline <- list(raw.conc = conc, raw.test=test, drID=drID, fit.conc=ytest$x, fit.test=ytest$y, spline = spltest, parameters=list("EC50"=EC.test[1], "yEC50"=yEC.test, "EC50.orig"=EC.orig[1], "yEC50.orig"=EC.orig[2]), fitFlag = fitFlag, reliable=NULL, control=control)

class(drFitSpline)<-"drFitSpline"

drFitSpline
}

