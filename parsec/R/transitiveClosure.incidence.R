transitiveClosure.incidence <-
function(m) {
    g <- incidence2cover(m)
    ct <- transitiveClosure.default(g)
    cover2incidence(ct)
}
