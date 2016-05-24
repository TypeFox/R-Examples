CADFtest <- function(model, X=NULL, type=c("trend", "drift", "none"), 
                     data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname=NULL, 
                     criterion=c("none", "BIC", "AIC", "HQC", "MAIC"), ...)
UseMethod("CADFtest")
 