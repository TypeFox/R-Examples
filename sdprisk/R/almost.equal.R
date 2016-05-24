almost.equal <- function(x, value) {
    vapply(Map(all.equal, x, rep.int(value, length(x))), isTRUE, logical(1L))
}
