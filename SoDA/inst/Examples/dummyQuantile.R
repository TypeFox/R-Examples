doQuantile <- function(oldQ, newdata, probs) {
    if(!exists("hiddenData")) {
        message("using dummy doQuantile--don't take the results seriously")
        hiddenData <<- numeric()
    }
    hiddenData <<- c(hiddenData, newdata)
    quantile(hiddenData, probs)
}
