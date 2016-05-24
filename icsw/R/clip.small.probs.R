clip.small.probs <- function(x, min.prob = NULL) {
    if (!is.null(min.prob)) {
        if(min(x) < min.prob) {
            warning("Small probabilities encountered. Replacing with min.prob.")
            x[x < min.prob] <- min.prob
        }
    }
    if (min(x) <= 0 & is.null(min.prob)) {
        warning("Probabilities <= 0 encountered and min.prob not set. Replacing with smallest value > 0.")
        x[x <= 0] <- min(x[x > 0])
    }
    return(x)
}