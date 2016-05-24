GHQ <- function(n, ndim, pruning=TRUE) {
#  list.of.packages <- c("statmod")
#  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#  if(length(new.packages)) install.packages(new.packages)
#  install.packages("statmod", dependencies=TRUE)
#if ( ! is.element('lme4', installed.packages()[,1]) ) install.packages('statmod')
#require(statmod)
quad <- statmod::gauss.quad(n=n, kind="hermite")
#quad <- gauss.quad(n=n, kind="hermite")
gride.nodes <- gride.weights <- list()
for (i in 1:ndim) {
  gride.nodes[[i]]   <- quad$nodes
  gride.weights[[i]] <- quad$weights }
Nodes   <- expand.grid(gride.nodes)
Weights <- expand.grid(gride.weights)
product <- apply(Weights,1,prod)
theta_m <- (quad$weights[1] * quad$weights[(n+1)/2]) / n^(ndim-1)
ind <- product > theta_m
if (pruning) {
Nodes <- Nodes[ind,]
Weights <- Weights[ind,]
product <- product[ind]
}
result <- list(nodes=Nodes, weights=Weights, product=product)
if (ndim == 1) lapply(result, as.matrix)
else result }
