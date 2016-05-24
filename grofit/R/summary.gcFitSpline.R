summary.gcFitSpline <-
function(object,...)
{

# object of class gcFitSpline

contents.fitted.spline  <- c("mu.spline", "lambda.spline", "A.spline", "integral.spline")

if ((is.na(object$fitFlag)==TRUE)|(object$fitFlag==FALSE)){
   table<-rep(NA,length(contents.fitted.spline))
}
else{
    table <- c(object$parameters$mu, object$parameters$lambda,  object$parameters$A, object$parameters$integral)
}

table               <- data.frame(t(table))
colnames(table)     <- contents.fitted.spline
summary.gcFitSpline <- table

}

