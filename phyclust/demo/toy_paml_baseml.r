library(phyclust, quiet = TRUE)

X <- seq.data.toy$org

### Fit a EE, JC69 model using emEM
EMC <- .EMControl(init.procedure = "emEM")
set.seed(1234)
K <- 4
ret.K <- phyclust(X, K, EMC = EMC)

### Obtain a ML tree for central sequences.
(ret.Mu <- paml.baseml(ret.K$Mu, opts = paml.baseml.control()))

### Construct an ancestral tree.
tree.anc <- read.tree(text = ret.Mu$best.tree)
plot(tree.anc, main = "tree of centers")

### Construct descent trees.
n.class <- ret.K$n.class
Tt <- ret.K$QA$Tt
tree.est <- tree.anc
for(k in 1:K){
  tree.dec <- gen.star.tree(n.class[k], total.height = Tt) 
  tree.est <- bind.tree(tree.est, tree.dec,
                          where = which(tree.est$tip.label == as.character(k)))
}
est.class <- rep(1:K, n.class)
plotnj(tree.est, X.class = est.class, main = "tree of sequences")

