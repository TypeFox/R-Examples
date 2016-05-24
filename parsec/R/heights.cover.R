heights.cover <-
function(z) {
    z <- cover2incidence(z)
    colSums(as.matrix(z))
}
