
##
require("slam")

##
x <- matrix(c(1L, 0L, 3L, 0L, 5L, 0L), ncol = 2,
    dimnames = list(1:3, LETTERS[1:2]))
x
s <- as.simple_triplet_matrix(x)

colapply_simple_triplet_matrix(s, identity)
rowapply_simple_triplet_matrix(s, identity)

s$v <- as.numeric(s$v)
simplify2array(colapply_simple_triplet_matrix(s, identity))

s$v <- as.complex(s$v)
simplify2array(colapply_simple_triplet_matrix(s, identity))

s$v <- as.list(s$v)
simplify2array(colapply_simple_triplet_matrix(s, identity))

s$v <- as.character(s$v)
simplify2array(colapply_simple_triplet_matrix(s, identity))

##
