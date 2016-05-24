is.preorder <-
function(m) all(binary(m), reflexivity(m), antisymmetry(m))
