weightedMean <- function(x, w) {
    .Call('IBHM_weightedMean', PACKAGE = 'IBHM', x, w)
}

weightedVar <- function(x, w) {
    .Call('IBHM_weightedVar', PACKAGE = 'IBHM', x, w)
}

weightedCov <- function(x, z, w) {
    .Call('IBHM_weightedCov', PACKAGE = 'IBHM', x, z, w)
}

weightedR <- function(y1, y2, w) {
    .Call('IBHM_weightedR', PACKAGE = 'IBHM', y1, y2, w)
}

tiedRanks <- function(x) {
    .Call('IBHM_tiedRanks', PACKAGE = 'IBHM', x)
}

weightedRho <- function(y1, y2, w) {
    .Call('IBHM_weightedRho', PACKAGE = 'IBHM', y1, y2, w)
}

