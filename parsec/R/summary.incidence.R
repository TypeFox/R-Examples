summary.incidence <-
function(object, ...) {
    number.of.elements <- nrow(object)
    number.of.comparability <- sum(object)*2 - number.of.elements
    list(number.of.elements=number.of.elements,
    number.of.comparability=number.of.comparability)
}
