.singleRaisedcos <- function(x, lo, center, hi) {
    if (x < lo || x > hi) {
        return(0)
    } else if (x < center) {
        if (lo == -Inf) {
            return(1)
        }
        return((cos((x - center) * pi / (center - lo)) + 1) / 2)
    } else if (x == center) {
        return(1)
    } else {
        if (hi == Inf) {
            return(1)
        }
        return((cos((x - center) * pi / (hi - center)) + 1) / 2)
    }
}


raisedcos <- function(x, lo, center, hi) {
    if (lo > center || lo > hi) {
        stop('"lo" must be the lower-bound of the interval <lo, hi>')
    }
    if (hi < center || hi < lo) {
        stop('"hi" must be the upper-bound of the interval <lo, hi>')
    }
    if (center < lo || center > hi) {
        stop('"center" must be within the interval <lo, hi>')
    }

    return(sapply(x, .singleRaisedcos, lo, center, hi))
}
