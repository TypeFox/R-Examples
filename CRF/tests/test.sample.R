library(CRF)

test.sample <- function(name, sample.method, dataset, cutoff=0.01, size=10000, ...)
{
  crf <- dataset$crf
  answer <- dataset$answer

  cat("  ", name, ": Sampling ... ", sep="")
  samples <- sample.method(crf, size, ...)
  
  samples.node.bel <- array(0, dim=c(crf$n.nodes, crf$max.state))
  for (i in 1:crf$n.nodes)
    for (j in 1:crf$max.state)
      samples.node.bel[i,j] <- sum(samples[,i] == j)
  samples.node.bel = samples.node.bel / rowSums(samples.node.bel)
  if (mean(abs(samples.node.bel - answer$node.bel)) < cutoff) {
    cat("Passed.\n")
  } else {
    cat("Failed ***\n")
    warning(name, ": Sampling may be incorrect!")
  }
}

cat("Testing dataset Small ...\n")
data(Small)
test.sample("Exact", sample.exact, Small)
test.sample("Chain", sample.chain, Small)
test.sample("Tree", sample.tree, Small)
test.sample("Cutset", sample.cutset, Small, 0.01, 10000, c(3))
test.sample("Cutset (Chain)", sample.cutset, Small, 0.01, 10000, c(1), "chain")
test.sample("JunctionTree", sample.junction, Small)
test.sample("Gibbs", sample.gibbs, Small, 0.01, 1000000, 10000)

cat("Testing dataset Chain ...\n")
data(Chain)
test.sample("Chain", sample.chain, Chain)
test.sample("Tree", sample.tree, Chain)
test.sample("Cutset", sample.cutset, Chain, 0.01, 10000, c(3))
test.sample("Cutset (Chain)", sample.cutset, Chain, 0.01, 10000, c(1), "chain")
test.sample("JunctionTree", sample.junction, Chain)
test.sample("Gibbs", sample.gibbs, Chain, 0.01, 1000000, 10000)

cat("Testing dataset Tree ...\n")
data(Tree)
test.sample("Tree", sample.tree, Tree)
test.sample("Cutset", sample.cutset, Tree, 0.01, 10000, c(3))
test.sample("JunctionTree", sample.junction, Tree)
test.sample("Gibbs", sample.gibbs, Tree, 0.01, 1000000, 10000)

cat("Testing dataset Loop ...\n")
data(Loop)
test.sample("Cutset", sample.cutset, Loop, 0.01, 10000, c(1,3))
test.sample("Cutset (Chain)", sample.cutset, Loop, 0.01, 10000, c(1), "chain")
test.sample("JunctionTree", sample.junction, Loop)
test.sample("Gibbs", sample.gibbs, Loop, 0.01, 1000000, 10000)

cat("Testing dataset Clique ...\n")
data(Clique)
test.sample("Cutset", sample.cutset, Clique, 0.01, 10000, c(1,3))
test.sample("JunctionTree", sample.junction, Clique)
test.sample("Gibbs", sample.gibbs, Clique, 0.01, 1000000, 10000)
