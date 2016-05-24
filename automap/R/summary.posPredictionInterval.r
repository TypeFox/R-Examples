summary.posPredictionInterval = function(object, ...)
{
    cat("Upper and lower boundary and the position relative to value:")
    print(summary(object$pos_prediction), ...)
    cat("\n")
    print(object["p"], ...)
    print(object["value"], ...)
}
