drBootSpline <-
function (conc, test, drID="undefined", control=grofit.control())
{

# /// conversion to handle even data.frame inputs
test <- as.vector(as.numeric(as.matrix(test)))
conc <- as.vector(as.numeric(as.matrix(conc)))

# /// check input values
if (is.vector(conc)==FALSE || is.vector(test)==FALSE) stop("Need concentration and treatment !")
if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
if (control$nboot.dr==0) stop("Number of bootstrap samples is zero! See grofit.control()")

if (control$neg.nan.act==FALSE){
   # /// removing missing data
   missings<-is.na(conc)|is.na(test)|!is.numeric(conc)|!is.numeric(test)
   conc<-conc[!missings]
   test<-test[!missings]

   # /// removing values wrongly assigned negative values
   negs<-(conc<0)|(test<0)
   conc<-conc[!negs]
   test<-test[!negs]
}
else{
   if (sum(is.na(conc)|is.na(test))) stop("NA values encountered. Program terminated")
   if ((sum((conc<0))>0)|(sum((test<0))>0)) stop("drFitSpline: Negative values encountered. Program terminated")
   if ((FALSE%in%is.numeric(conc))||(FALSE%in%is.numeric(test))) stop("drFitSpline: Non numeric values encountered. Program terminated")
}

# /// check length of data after removal of values
if (length(test)<6){
	warning("drBootSpline: There is not enough valid data. Must have at least 6 unique values!")
        drBootSpline <- list(raw.conc=conc, raw.test=test, drID=drID, boot.conc=NA, boot.test=NA, boot.drSpline = NA, ec50.boot=NA, bootFlag=FALSE, control=control)
        class(drBootSpline) <- "drBootSpline"
        return(drBootSpline)
}

if (length(test)<control$have.atleast){
	warning("drBootSpline: number of valid data points is below the number specified in 'have.atleast'. See grofit.control().")
        drBootSpline <- list(raw.conc=conc, raw.test=test, drID=drID, boot.conc=NA, boot.test=NA, boot.drSpline = NA, ec50.boot=NA, bootFlag=FALSE, control=control)
        class(drBootSpline) <- "drBootSpline"
        return(drBootSpline)
}


# /// transformation of data...
if (control$log.x.dr==TRUE) conc.log <- log(conc+1)
if (control$log.y.dr==TRUE) test.log <- log(test+1)
if (control$log.x.dr==TRUE) conc.boot<- log(conc+1) else conc.boot<-conc
if (control$log.y.dr==TRUE) test.boot<- log(test+1) else test.boot<-test

# /// Initialize some variables
boot.x      <- array(NA,c(control$nboot.dr,1000))
boot.y      <- array(NA,c(control$nboot.dr,1000))
ECtest.boot <- seq(0,0,length.out=control$nboot.dr)
splinefit   <- list()
sa          <- seq(1,length(conc.boot))

# /// begin bootstrapping
for(b in 1:control$nboot.dr){
	s      <- sample(sa, length(conc.boot), replace = TRUE)
	s.conc <- conc.boot[s]
	# ensure to choose enough values for spline smoothing
	while (length(unique(s.conc))<5)
	{
            s      <- sample(sa, length(conc.boot), replace = TRUE)
            s.conc <- conc.boot[s]
      	}
        s.test                            <- test.boot[s];
	spltest                           <- NULL
        control.changed                   <- control
        control.changed$suppress.messages <- TRUE

	# /// call spline fit (suppressed messages)
        splinefit[[b]] <- drFitSpline(s.conc, s.test, drID, control.changed)
        spltest        <- splinefit[[b]]$spline

	# /// create tables of fitted values from spline fit
        boot.x[b,1:length(splinefit[[b]]$fit.conc)] <- splinefit[[b]]$fit.conc
        boot.y[b,1:length(splinefit[[b]]$fit.test)] <- splinefit[[b]]$fit.test

	if (is.null(spltest)==TRUE){
            cat("Spline could not be fitted to data!\n")
		if (is.null(control$smooth.dr)==TRUE){
	 	    cat("This might be caused by usage of smoothing parameter NULL\n")
		}
	    stop("Error in drBootSpline")	
	}
        ECtest.boot[b] <- splinefit[[b]]$parameters$EC50;
} # of for(b in 1:control$nboot.dr){

# /// remove infinite values (which probably will not occure)
ECtest.boot[which(!is.finite(ECtest.boot))] <- NA

# /// remove negative values
if (control$clean.bootstrap == TRUE) ECtest.boot[which(ECtest.boot<0)] <- NA

# /// mean and st.dev. of bootstrap sample
m.test <- mean(ECtest.boot, na.rm=TRUE)
s.test <- sd(ECtest.boot, na.rm=TRUE)

if (control$suppress.messages==FALSE){
	cat("=== Bootstrapping of dose response curve ==========\n")
	cat("--- EC 50 -----------------------------------------\n")
	cat("\n")
	cat(paste("Mean  : ", as.character(m.test), "StDev : ", as.character(s.test),"\n"))
	cat(paste("90% CI: ", as.character(c(m.test-1.645*s.test,m.test+1.645*s.test))))
	cat("\n")
	cat(paste("95% CI: ", as.character(c(m.test-1.96*s.test,m.test+1.96*s.test)))) 
	cat("\n\n")
}

EC50 <- data.frame("EC50.boot"=m.test, "EC50.sd"=s.test, "CI.90.lo"=m.test-1.645*s.test, "CI.90.up"=m.test+1.645*s.test, "CI.95.lo"=m.test-1.96*s.test, "CI.95.up"=m.test+1.96*s.test)

if (control$log.x.dr==TRUE && control$suppress.messages==FALSE){
	cat("\n")
	cat("--- EC 50 in original scale -----------------------\n")
	cat("\n")
	cat(paste("Mean  : ", as.character(exp(m.test)-1),"\n"))
	cat(paste("90% CI: ", as.character(c(exp(m.test-1.645*s.test)-1,exp(m.test+1.645*s.test)-1))))
	cat("\n")
	cat(paste("95% CI: ", as.character(c(exp(m.test-1.96*s.test)-1,exp(m.test+1.96*s.test)-1))))
	cat("\n\n")
}

drBootSpline        <- list(raw.conc=conc, raw.test=test, drID=drID, boot.conc=boot.x, boot.test=boot.y, boot.drSpline = splinefit, ec50.boot=ECtest.boot, bootFlag=TRUE, control=control)
class(drBootSpline) <- "drBootSpline"
drBootSpline

}

