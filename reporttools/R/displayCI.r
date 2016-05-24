displayCI <- function(ci, digit = 2, unit = "", text = "none"){
    ci <- format(round(ci, digit), nsmall = digit)
    d <- 1
    if (is.matrix(ci) == TRUE){d <- nrow(ci)} else {ci <- matrix(ci, ncol = 2)}
    ci <- sub(' ', '', ci)    # eliminate spaces
    disp.ci <- rep(NA, d)
    for (i in 1:d){
        if (identical(text, "none")){string <- paste("[", ci[i, 1], unit, ", ", ci[i, 2], unit, "]", sep = "")}
        if (identical(text, "german")){string <- paste("von ", ci[i, 1], unit, " bis ", ci[i, 2], unit, sep = "")}
        if (identical(text, "english")){string <- paste("from ", ci[i, 1], unit, " to ", ci[i, 2], unit, sep = "")}
        disp.ci[i] <- string
        }
    return(disp.ci)
}
