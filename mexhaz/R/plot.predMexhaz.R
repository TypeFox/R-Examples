plot.predMexhaz <- function(x,which=c("surv","hazard"),conf.int=TRUE,...){
    which <- match.arg(which)
    if (x$type=="multiobs"){
        stop("The 'plot' function applies only to predictions realised on a single vector of covariables.")
    }
    time.pts <- x$results$time.pts
    if (which=="hazard"){
        plot(time.pts,x$results$hazard,type="l",xaxs="i",xlab="Time",ylab="Hazard",...)
    }
    else {
        plot(c(0,time.pts),c(1,x$results$surv),type="l",xaxs="i",xlab="Time",ylab="Survival",...)
    }
    if (conf.int==TRUE & x$ci.method%in%c("delta","simul")){
        if (which=="hazard"){
            points(time.pts,x$results$hazard.inf,type="l",lty="dashed",...)
            points(time.pts,x$results$hazard.sup,type="l",lty="dashed",...)
        }
        else {
            points(c(0,time.pts),c(1,x$results$surv.inf),type="l",lty="dashed",...)
            points(c(0,time.pts),c(1,x$results$surv.sup),type="l",lty="dashed",...)
        }
    }
    if (conf.int==TRUE & !(x$ci.method%in%c("delta","simul"))){
        warning("Confidence limits could not be displayed because they were missing from the predMexhaz object.")
    }
}
