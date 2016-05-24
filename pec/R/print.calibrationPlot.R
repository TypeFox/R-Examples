### print.calibrationPlot.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct  4 2015 (09:49) 
## Version: 
## last-updated: Oct  5 2015 (10:32) 
##           By: Thomas Alexander Gerds
##     Update #: 28
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @export
print.calibrationPlot <- function(x,...){
    if (x$model.type=="survival"){
        estring <- paste0("\n - a total of ",
                          x$summary$Event[1],
                          " were observed to have an event,\n - a total of ",
                          x$summary$Lost,
                          " were lost to follow-up.")
    }else{
         estring <- paste0("\n - a total of ",
                           x$summary$Event[1],
                           " were observed to have the event of interest (cause: ",x$cause,"),\n - a total of ",
                           sum(x$summary$Event[-1]),
                           " had a competing risk ",
                           " \n - a total of ",
                           x$summary$Lost,
                           " were lost to follow-up.")
     }
    cat("\nCalibration of ",ifelse(x$model.type=="survival" && x$type=="survival","survival","risk")," predictions for ",
        x$summary$n,
        " subjects.\n\nUntil time ", x$time,
        " a total of ",x$summary$Event.free,
        " were observed event-free, ",
        estring,
        "\n",
        sep="")
    if (x$method=="quantile"){
        cat("\nAverage predictions and outcome in prediction quantiles:\n\n")
        print(x$plotFrames,...)
    }else{
         cat("\nSummary of predictions and outcome:\n")
         nix <- lapply(1:length(x$plotFrames),
                       function(n){
                           cat("\n",names(x$plotFrames)[[n]],":\nOutcome:\n",sep="")
                           print(summary(x$plotFrames[[1]][,1]))
                           cat("\nPredictions:\n",sep="")                           
                           print(summary(x$plotFrames[[1]][,2]))
                       })
     }
    cat("\nOutcome frequencies (Obs) were obtained with the ",
        ifelse(x$pseudo,
               "jackknife pseudo-value method ",
               ifelse(x$model.type=="survival","Kaplan-Meier method.","Aalen-Johansen method.")),
        "\n",
        sep="")
}



#----------------------------------------------------------------------
### print.calibrationPlot.R ends here
