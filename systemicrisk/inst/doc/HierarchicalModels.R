## ------------------------------------------------------------------------
library(systemicrisk)

## ------------------------------------------------------------------------
l <- c(714,745,246, 51,847)
a <- c(872, 412, 65, 46,1208)

## ------------------------------------------------------------------------
mod <- Model.additivelink.exponential.fitness(n=5,alpha=-2.5,beta=0.3,gamma=1.0,
                 lambdaprior=Model.fitness.genlambdaparprior(ratescale=500))

## ------------------------------------------------------------------------
thin <- choosethin(l=l,a=a,model=mod,silent=TRUE)
thin

## ------------------------------------------------------------------------
res <- sample_HierarchicalModel(l=l,a=a,model=mod,nsamples=1e3,thin=thin,silent=TRUE)

## ------------------------------------------------------------------------
res$L[[1]]
res$L[[2]]

## ----fig.width=7,fig.height=4--------------------------------------------
plot(ecdf(sapply(res$L,function(x)x[1,2])))

## ------------------------------------------------------------------------
diagnose(res)

## ----fig.width=7,fig.height=4--------------------------------------------
plot(sapply(res$L,function(x)x[1,2]),type="b")

## ----fig.width=7,fig.height=4--------------------------------------------
plot(res$theta[1,],type="b")

## ----fig.width=7,fig.height=4--------------------------------------------
acf(sapply(res$L,function(x)x[1,2]))

