is.partialorder <-
function(m) all(is.preorder(m), transitivity(m))
