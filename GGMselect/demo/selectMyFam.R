p=30
n=30
# generate graph
eta=0.11
Gr <- simulateGraph(p,eta)
# generate data
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)
# generate a family of candidate graphs with glasso
library("glasso")
MyFamily <- NULL
for (j in 1:3){
  MyFamily[[j]] <- abs(sign(glasso(cov(X),rho=j/5)$wi))
  diag(MyFamily[[j]]) <- 0
}
# select a graph within MyFamily
GMF <- selectMyFam(X,MyFamily)
# plot the result
library(network)
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
gV <- network(Gr$G)
a<-plot(gV, usearrows = FALSE)
gMyFam <- network(GMF$G)
plot(gMyFam, coord=a, usearrows = FALSE)
par(old.par)
