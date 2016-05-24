
library("slam")
set.seed(201311)

###
x <- matrix(rnorm(100), nrow = 20,
    dimnames =  list(1:20, LETTERS[1:5])
)
x[sample(100, 80)] <- 0

s <- as.simple_triplet_matrix(x)
s

##
identical(apply(x, 2L, var), colapply_simple_triplet_matrix(s, var))
identical(apply(x, 1L, var), rowapply_simple_triplet_matrix(s, var))

##
k <- 1:2
z <- var(x[, k], x[, -k])
identical(z, crossapply_simple_triplet_matrix(s[, k], s[, -k], FUN = var))
identical(z, crossapply_simple_triplet_matrix(x[, k], s[, -k], FUN = var))

identical(z, 
    tcrossapply_simple_triplet_matrix(t(s[, k]), t(s[, -k]), FUN = var))
identical(z, 
    tcrossapply_simple_triplet_matrix(t(x[, k]), t(s[, -k]), FUN = var))

z <- var(x)
identical(z, crossapply_simple_triplet_matrix(s, FUN = var))

## null-dimensions
z <- var(x[, 0], x)
z
all.equal(z, crossapply_simple_triplet_matrix(s[, 0], s, FUN = var))
all.equal(z, crossapply_simple_triplet_matrix(x[, 0], s, FUN = var))

try(crossapply_simple_triplet_matrix(x[, 0], s, FUN = var, use = "all.obs"))

z <- var(x, x[, 0])
z
all.equal(z, crossapply_simple_triplet_matrix(s, s[, 0], FUN = var))
all.equal(z, crossapply_simple_triplet_matrix(x, s[, 0], FUN = var))


z <- var(x[, 0])
z
all.equal(z, crossapply_simple_triplet_matrix(s[, 0], s[, 0], FUN = var))
all.equal(z, crossapply_simple_triplet_matrix(x[, 0], s[, 0], FUN = var))

all.equal(z, crossapply_simple_triplet_matrix(s[, 0], FUN = var))

z <- var(x[0, ])
z
all.equal(z, crossapply_simple_triplet_matrix(s[0, ], s[0, ], FUN = var))
all.equal(z, crossapply_simple_triplet_matrix(x[0, ], s[0, ], FUN = var))

all.equal(z, crossapply_simple_triplet_matrix(s[0, ], FUN = var))

## non-scalar
z <- crossapply_simple_triplet_matrix(s, s, FUN = ">")
all.equal(z, crossapply_simple_triplet_matrix(x, s, FUN = ">"))

all.equal(z[lower.tri(z)],
    crossapply_simple_triplet_matrix(s, FUN = ">")[lower.tri(z)])

###
