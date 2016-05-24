library(phyclust, quiet = TRUE)

### Examples to use ms().
set.seed(1234)
K <- 4                                # Number of clusters.
N <- 50                               # Number of sequences.
Eta <- 0.05 + runif(K)
Eta <- Eta / sum(Eta)                 # Population proportions.
rate.anc <- 0.1                       # r_a.
rate.dec <- 0.1                       # r_d.

### Generate an ancestral tree.
ms.anc <- ms(K, opts = paste("-T -G", rate.anc, sep = " "))
tree.anc <- read.tree(text = ms.anc[3])
tree.anc$tip.label <- paste("a", 1:K, sep = "")

### Generate descent trees and attach them to the ancestral tree.
N.k <- as.vector(rmultinom(1, N, Eta))
ms.dec <- NULL
tree.dec <- NULL
tree.joint <- tree.anc
for(k in 1:K){
  ms.dec[[k]] <- ms(N.k[k], opts = paste("-T -G", rate.dec, sep = " "))
  tree.dec[[k]] <- read.tree(text = ms.dec[[k]][3])
  tree.dec[[k]]$tip.label <- paste("d", k, ".", 1:N.k[k], sep = "")
  tree.joint <- bind.tree(tree.joint, tree.dec[[k]],
                          where = which(tree.joint$tip.label ==
                                        paste("a", k, sep = "")))
}

### Plot the trees.
par(mfrow = c(2, 3))
plot(tree.anc, main = paste("anc (", K, ")", sep = ""))
for(k in 1:K){
  plot(tree.dec[[k]], main = paste("dec", k, " (", N.k[k], ")", sep = ""))
}
plot(tree.joint, main = paste("joint (", N, ")", sep = ""))
