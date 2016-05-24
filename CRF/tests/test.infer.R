library(CRF)

test.infer <- function(name, infer.method, dataset, cutoff=1e-8, ...)
{
  crf <- dataset$crf
  answer <- dataset$answer

  cat("  ", name, ": Inferring ... ", sep="")
  belief <- infer.method(crf, ...)
  
  node.error <- sapply(1:crf$n.nodes, function(i) max(abs(belief$node.bel[i,] - answer$node.bel[i,])))
  edge.error <- sapply(1:crf$n.edges, function(i) max(abs(belief$edge.bel[[i]] - answer$edge.bel[[i]])))
  if (max(abs(c(node.error, edge.error, belief$logZ - answer$logZ))) < cutoff) {
    cat("Passed.\n")
  } else {
    cat("Failed ***\n")
    warning(name, ": Inference is incorrect!")
  }
}

cat("Testing dataset Small ...\n")
data(Small)
test.infer("Exact", infer.exact, Small)
test.infer("Chain", infer.chain, Small)
test.infer("Tree", infer.tree, Small)
test.infer("Cutset", infer.cutset, Small, 1e-8, c(3))
test.infer("Cutset (Chain)", infer.cutset, Small, 1e-8, c(1), "chain")
test.infer("JunctionTree", infer.junction, Small)
test.infer("Sample", infer.sample, Small, 1e-2, sample.exact, 1000000)
test.infer("LBP", infer.lbp, Small)
test.infer("TRBP", infer.trbp, Small)

cat("Testing dataset Chain ...\n")
data(Chain)
test.infer("Chain", infer.chain, Chain)
test.infer("Tree", infer.tree, Chain)
test.infer("Cutset", infer.cutset, Chain, 1e-8, c(3))
test.infer("Cutset (Chain)", infer.cutset, Chain, 1e-8, c(1), "chain")
test.infer("JunctionTree", infer.junction, Chain)
test.infer("Sample", infer.sample, Chain, 1e-2, sample.chain, 1000000)
test.infer("LBP", infer.lbp, Chain)
test.infer("TRBP", infer.trbp, Chain)

cat("Testing dataset Tree ...\n")
data(Tree)
test.infer("Tree", infer.tree, Tree)
test.infer("Cutset", infer.cutset, Tree, 1e-8, c(3))
test.infer("JunctionTree", infer.junction, Tree)
test.infer("Sample", infer.sample, Tree, 1e-2, sample.tree, 1000000)
test.infer("LBP", infer.lbp, Tree)
test.infer("TRBP", infer.trbp, Tree)

cat("Testing dataset Loop ...\n")
data(Loop)
test.infer("Cutset", infer.cutset, Loop, 1e-8, c(3))
test.infer("Cutset (Chain)", infer.cutset, Loop, 1e-8, c(1), "chain")
test.infer("JunctionTree", infer.junction, Loop, 0.01)
test.infer("Sample", infer.sample, Loop, 1e-2, sample.exact, 1000000)
test.infer("LBP", infer.lbp, Loop, 0.01)
test.infer("TRBP", infer.trbp, Loop, 0.1)

cat("Testing dataset Clique ...\n")
data(Clique)
test.infer("Cutset", infer.cutset, Clique, 1e-8, c(1,3))
test.infer("JunctionTree", infer.junction, Clique, 0.1)
test.infer("Sample", infer.sample, Clique, 1e-2, sample.exact, 1000000)
test.infer("LBP", infer.lbp, Clique, 0.1)
test.infer("TRBP", infer.trbp, Clique, 0.1)
