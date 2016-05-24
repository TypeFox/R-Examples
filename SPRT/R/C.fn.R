C.fn <-
function(distribution, h) {
    if (distribution == "normal") {
        h
    } else if (distribution == "bernoulli") {
        log(h / (1-h))
    } else if (distribution == "poisson") {
        log(h)
    } else if (distribution == "exponential") {
        -1/h
    }
}
