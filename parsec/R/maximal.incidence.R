maximal.incidence <-
function(z) {
    tot <- rowSums(as.matrix(z))
    tot == 1
}
