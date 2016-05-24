
twoTimes <- function(x) {
    .Call(twoTimesImpl, x)   # twoTimesImpl is declared in NAMESPACE, and defined in src/
}

