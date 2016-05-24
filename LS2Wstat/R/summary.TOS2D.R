summary.TOS2D <-function (object, ...) {


    cat("\n2D bootstrap test of stationarity\n") 
    cat("      object of class TOS2D\n----------------------------------\n")
    cat("\nsummary\n=======\n")
    cat("data:", object$data.name,"\n")
    cat("Observed test statistic:",round(object$samples[1],3),"\n")
    cat("bootstrap p-value:",round(object$p.value,3),"\n\n") 

}
