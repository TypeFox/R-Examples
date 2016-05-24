library(phyclust, quiet = TRUE)

### Examples to apply ms() and seqgen() for nucleotide data.
set.seed(1234)
K <- 4                                # Number of clusters.
N <- 200                              # Number of sequences.
L <- 400                              # Number of sites.
Eta <- 0.05 + runif(K)
Eta <- Eta / sum(Eta)                 # Population proportions.
pi <- c(0.35, 0.17, 0.21, 0.27)       # Equilibrium probabilities.
kappa <- 9.0                          # For HKY85.
rate.anc <- 0.1                       # r_a.
rate.dec <- 0.1                       # r_d.
rate.scale <- 0.10                    # h_s.

### Generate trees by ms().
N.K <- as.vector(rmultinom(1, N - K, Eta)) + 1
tree.K <- gen.unit.K(K, N.K, rate.anc, rate.dec)

### Generate sequences by seqgen().
seq.K.equal <- gen.seq.HKY(tree.K$equal, pi, kappa, L, rate.scale)
da.K.equal <- read.seqgen(seq.K.equal)

### Generate sequences givn ancestor sequences by seqgen().
rooted.tree <- tree.K$anc
anc.seq <- da.K.equal$org[1,]
new <- gen.seq.HKY(rooted.tree, pi, kappa, L, rate.scale, anc.seq = anc.seq)
