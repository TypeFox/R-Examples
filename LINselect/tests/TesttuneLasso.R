#source("charge.R")
library("LINselect")

set.seed(15) # pour pouvoir comparer des exécutions sucessives
ex <- simulData(p=100,n=50,r=0.95,rSN=10)

res.tuneLasso <- tuneLasso(ex$Y,ex$X,method="Glasso",verbose=TRUE)
res.tuneLasso$Glasso$CV$support
res.tuneLasso$Glasso$Ls$support

res.tuneLasso <- tuneLasso(ex$Y,ex$X,method="lasso",verbose=TRUE)
res.tuneLasso$lasso$CV$support
res.tuneLasso$lasso$Ls$support

## Time consuming
# set.seed(5) # pour pouvoir comparer des exécutions sucessives
# ex <- simulData(p=100,n=100,r=0.95,rSN=10)
# res.tuneLasso <- tuneLasso(ex$Y,ex$X, dmax=20)
# 
# 
# set.seed(67) # pour pouvoir comparer des exécutions sucessives
# ex <- simulData(p=50,n=1000,r=0.95,rSN=10)
# res.tuneLasso <- tuneLasso(ex$Y,ex$X, dmax=20)
# res.tuneLasso$lasso$CV$coeff
# res.tuneLasso$lasso$Ls$coeff
# res.tuneLasso$Glasso$CV$support
# res.tuneLasso$Glasso$Ls$support



