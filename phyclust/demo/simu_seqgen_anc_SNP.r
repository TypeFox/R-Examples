library(phyclust, quiet = TRUE)

### Examples to apply ms() and seq.gen() for SNP data.
set.seed(1234)
N <- 1001                             # Number of sequences.
L <- 8                                # Number of sites.
pi <- c(0.5, 0.5)                     # Equilibrium probabilities.
Tt <- 0.48                            # Evolved time.

### Generate a star tree with N tips and followed by a sequences generation.
star.tree <- gen.star.tree(N)
seq.star <- gen.seq.SNP(star.tree, pi, L, Tt)
