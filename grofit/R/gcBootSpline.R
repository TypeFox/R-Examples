gcBootSpline <-
function (time, data, gcID ="undefined", control=grofit.control())
{

# /// check input parameters
if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
if (control$nboot.gc==0) stop("Number of bootstrap samples is zero! See grofit.control()")

# /// conversion to handle even data.frame inputs
time <- as.vector(as.numeric(as.matrix(time)))
data <- as.vector(as.numeric(as.matrix(data)))

# /// check length of input data
if (length(time)!=length(data)) stop("gcBootSpline: length of input vectors differ!")

# /// determine which values are not valid
bad.values <- (is.na(time))|(time<0)|(is.na(data))|(data<0)|(!is.numeric(time))|(!is.numeric(data))

# /// remove bad values or stop program
if (TRUE%in%bad.values){
    if (control$neg.nan.act==FALSE){
        time    <- time[!bad.values]
        data    <- data[!bad.values]
    }
    else stop("Bad values in gcBootSpline")
}

# /// check length of data after removal of values
if (length(data)<6){
	warning("gcBootSpline: There is not enough valid data. Must have at least 6 unique values!")
        gcBootSpline <- list(raw.time=time, raw.data=data, gcID =gcID, boot.time=NA, boot.data=NA, boot.gcSpline=NA, lambda=NA, mu=NA, A=NA, integral=NA, bootFlag=FALSE, control=control)
        class(gcBootSpline) <- "gcBootSpline"
        return(gcBootSpline)
}

# /// apply transformation
if (control$log.x.gc==TRUE) time <- log(1+time)
if (control$log.y.gc==TRUE) data <- log(1+data)

# /// initialize variables
A           <- NA
mu          <- NA
lambda      <- NA
integral    <- NA
boot.y      <- array(NA,c(control$nboot.gc,length(time)))
boot.x      <- array(NA,c(control$nboot.gc,length(time)))
nonpara     <- list()

# /// set fit option temporally to spline fit to pass the check
#     of the input parameters in gcFitSpline.R.
#     Otherwise setting of control$fit.opt="m" and control$nboot.gc>0
#     will cause an error in gcFitSpline.R.
control.change <- control
control.change$fit.opt <- "s"

if (control$nboot.gc>0){
	# /// Begin bootstrapping ---------------------------------------------------------------------------------------------
	for (j in 1:control$nboot.gc){
		choose <- sample(1:length(time),length(time),replace=TRUE)
		# /// make sure to choose enough values
		while (length(unique(choose))<5 ){
			choose <- sample(1:length(time),length(time),replace=TRUE)
		}

		time.cur <- time[choose]
		data.cur <- data[choose]
		
		nonpara[[j]] <- gcFitSpline(time.cur, data.cur, gcID, control.change)
		boot.y[j,1:length(nonpara[[j]]$fit.data)]  <- nonpara[[j]]$fit.data
		boot.x[j,1:length(nonpara[[j]]$fit.time)]  <- nonpara[[j]]$fit.time
		
		lambda[j]   <- nonpara[[j]]$parameters$lambda
		mu[j]       <- nonpara[[j]]$parameters$mu
		A[j]        <- nonpara[[j]]$parameters$A
		integral[j] <- nonpara[[j]]$parameters$integral
	
	} # /// end of bootstrap loop
	
	# /// set infinte values to NA
	lambda[which(!is.finite(lambda))]   <- NA
	mu[which(!is.finite(lambda))]       <- NA
	A[which(!is.finite(lambda))]        <- NA
	integral[which(!is.finite(lambda))] <- NA
	
	# /// remove bad entries in parameters
	if (control$clean.bootstrap==TRUE){
		lambda[which(lambda<0)]     <- NA
		mu[which(mu<0)]             <- NA
		A[which(A<0)]               <- NA
		integral[which(integral<0)] <- NA
	}

} # of if (control$nboot.gc>0){

gcBootSpline <- list(raw.time=time, raw.data=data, gcID =gcID, boot.time=boot.x, boot.data=boot.y, boot.gcSpline=nonpara, lambda=lambda, mu=mu, A=A, integral=integral, bootFlag=TRUE, control=control)

class(gcBootSpline) <- "gcBootSpline"
gcBootSpline
}

