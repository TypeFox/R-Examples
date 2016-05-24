## Test the hypothesesTree creation class.

set.seed(130229)
library('structSSI')
library('ape')
library('igraph')
tree.1 <- as.igraph(rtree(10))
tree.2 <- as.igraph(rtree(50))

V(tree.1)$name <- paste("hyp", c(1:19))
V(tree.2)$name <- paste("hyp", c(1:99))
tree.1.el <- get.edgelist(tree.1)
tree.2.el <- get.edgelist(tree.2)

unadjp.1 <- c(runif(5, 0, 0.01), runif(14, 0, 1))
names(unadjp.1) <- paste("hyp", c(1:19))
unadjp.2 <- c(runif(10, 0.01), runif(89, 0, 1))
names(unadjp.2) <- paste("hyp", c(1:99))

## The  hierarchical adjustment procedure
## applied to this class.
adjust1 <- hFDR.adjust(unadjp.1, tree.1.el)
adjust2 <- hFDR.adjust(unadjp.2, tree.2.el) # Correctly gives warning.

## Can plot results for the tree without the warning.
plot(adjust1)

