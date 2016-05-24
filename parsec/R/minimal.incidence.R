minimal.incidence <-
function(z) {
    tot <- colSums(as.matrix(z))
    tot == 1
}
