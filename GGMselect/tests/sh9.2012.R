library(GGMselect)

p=50
n=60
eta=0.092
extraeta=eta/5
Gr <- list(NULL)
iG=1
set.seed(iG)
Gr[[iG]] <- simulateGraph(p,eta,extraeta=extraeta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr[[iG]]$C)
X <- scale(X,center=TRUE,scale=FALSE)
#Gr.Chap10 <-  selectFast(X, verbose=T, dmax=10, K=2.5, family="C01")
#Gr.Chap25 <-  selectFast(X, verbose=T,dmax=25, K=2.5, family="C01")
print( selectFast(X,verbose=T, dmax=26, K=2.5, family="C01"))
# gc(verbose=FALSE)
# -----------------------------------------------------
p=50
n=1000
eta=0.092
extraeta=eta/5
Gr <- list(NULL)
iG=1
set.seed(iG)
Gr[[iG]] <- simulateGraph(p,eta,extraeta=extraeta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr[[iG]]$C)
X <- scale(X,center=TRUE,scale=FALSE)
print(  selectFast(X, dmax=10, K=2.5, family="C01"))
print( selectFast(X, dmax=25, K=2.5, family="C01"))
print( selectFast(X, dmax=26, K=2.5, family="C01"))
# -----------------------------------------------------
p=60
n=60
eta=0.092
extraeta=eta/5
Gr <- list(NULL)
iG=1
set.seed(iG)
Gr[[iG]] <- simulateGraph(p,eta,extraeta=extraeta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr[[iG]]$C)
X <- scale(X,center=TRUE,scale=FALSE)
print( selectFast(X,verbose=T,  dmax=10, K=2.5, family="C01"))
print(selectFast(X, verbose=T,dmax=25, K=2.5, family="C01"))
print(selectFast(X, dmax=26, K=2.5, family="C01"))
# -----------------------------------------------------
p=60
n=1000
eta=0.092
extraeta=eta/5
Gr <- list(NULL)
iG=1
set.seed(iG)
Gr[[iG]] <- simulateGraph(p,eta,extraeta=extraeta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr[[iG]]$C)
X <- scale(X,center=TRUE,scale=FALSE)
print(selectFast(X, dmax=10, K=2.5, family="C01"))
print(selectFast(X, dmax=25, K=2.5, family="C01"))
print(selectFast(X, dmax=26, K=2.5, family="C01"))
# -----------------------------------------------------
p=60
n=1000
eta=0.092
extraeta=eta/5
Gr <- list(NULL)
iG=1
set.seed(iG)
Gr[[iG]] <- simulateGraph(p,eta,extraeta=extraeta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr[[iG]]$C)
X <- scale(X,center=TRUE,scale=FALSE)
print(selectFast(X, dmax=26, K=2.5, family="C01"))

