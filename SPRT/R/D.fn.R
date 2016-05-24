D.fn <-
function(distribution, h) {
    if (distribution == "normal") {
        (h^2)/2
    } else if (distribution == "bernoulli") {
        -1 * log(1-h)
    } else if (distribution == "poisson") {
        h
    } else if (distribution == "exponential") {
        log(h)
    }
}
