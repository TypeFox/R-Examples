library(GGMselect)


p=50
n=2000

K <- c(1,1.25,1.5,1.75,2,2.5,5)
 dmax <- p-10

eta=0.092
extraeta=eta/5


Gr <- list(NULL)
iG=1
  set.seed(iG)
  Gr[[iG]] <- simulateGraph(p,eta,extraeta=extraeta)
iS=12
    set.seed(iS*(pi/3.1415)**iG)
    X <- rmvnorm(n, mean=rep(0,p), sigma=Gr[[iG]]$C)
#X <- scale(X,center=TRUE,scale=FALSE)

#calcModC01(X)
print( selectFast(X, dmax=dmax, K=K, family="C01") )

