.singleTriangle <- function(x, lo, center, hi) {
    if (x < center) {
        if (lo == -Inf) {
            return(1)
        }
        return(pmax(0, (x - lo) / (center - lo)))
    } else if (x == center) {
        return(1)
    } else {
        if (hi == Inf) {
            return(1)
        }
        return(pmax(0, (hi - x) / (hi - center)))
    }
}


triangle <- function(x, lo, center, hi) {
    if (lo > center || lo > hi) {
        stop('"lo" must be the lower-bound of the interval <lo, hi>')
    }
    if (hi < center || hi < lo) {
        stop('"hi" must be the upper-bound of the interval <lo, hi>')
    }
    if (center < lo || center > hi) {
        stop('"center" must be within the interval <lo, hi>')
    }

    return(sapply(x, .singleTriangle, lo, center, hi))
}
