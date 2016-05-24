p=30
n=30
# simulate graph
eta=0.11
Gr <- simulateGraph(p,eta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)
# estimate graph
GRest <- selectFast(X, family="C01")
# Neighb and G are 2 forms of the same result
a <- convertGraph(GRest$C01$Neighb)
cat("Is G equal to Neighb?\n")
print(all.equal(a, GRest$C01$G)) # TRUE
# recalculate the graph by MyFamily
GMF <- selectMyFam(X, list(a))
cat("Is G the same?\n")
print(all.equal(a,GMF$G)) # TRUE
