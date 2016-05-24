depths.cover <-
function(z) {
    z <- cover2incidence(z)
    rowSums(as.matrix(z))
}
