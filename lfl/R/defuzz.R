defuzz <- function(degrees,
                   values,
                   type=c('mom', 'fom', 'lom', 'dee')) {
    if (!is.vector(degrees) || !is.numeric(degrees)) {
        stop("'degrees' must be a numeric vector")
    }
    if (!is.vector(values) || !is.numeric(values)) {
        stop("'values' must be a numeric vector")
    }
    if (length(degrees) != length(values)) {
        stop("The length of 'degrees' and 'values' must be the same")
    }
    if (length(degrees) < 3) {
        stop("The length of 'degrees' must be at least 3")
    }

    if (min(degrees) < 0 || max(degrees) > 1) {
        stop("Values of 'degrees' must be truth values in the interval [0,1]")
    }

    type <- match.arg(type)
    i <- which(degrees == max(degrees))

    if (type == 'dee') {
        l <- length(degrees)
        diff <- degrees[2:l] - degrees[1:(l-1)]
        if (all(diff >= 0)) {
            type <- 'fom'
        } else if (all(diff <= 0)) {
            type <- 'lom'
        } else {
            type <- 'mom'
        }
    }

    if (type == 'fom') {
        return(values[min(i)])
    } else if (type == 'lom') {
        return(values[max(i)])
    } else { # if (type == 'mom') {
        return(mean(values[i]))
    }
}

