print.NHEMOtree <-
function(x, ...){
    cat("S metric:                    ", x$S_Metric, "\n")
    cat("Misclassification (min/max): ", x$Misclassification_total, "\n")
    cat("Costs (min/max):             ", x$Costs_total, "\n")
    invisible(x)
}
