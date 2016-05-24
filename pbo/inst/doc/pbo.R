## ------------------------------------------------------------------------
set.seed(765)
n <- 100
t <- 2400
m <- data.frame(matrix(rnorm(n*t),nrow=t,ncol=n,dimnames=list(1:t,1:n)),
                check.names=FALSE)

sr_base <- 0
mu_base <- sr_base/(252.0)
sigma_base <- 1.00/(252.0)**0.5

for ( i in 1:n ) {
  m[,i] = m[,i] * sigma_base / sd(m[,i]) # re-scale
  m[,i] = m[,i] + mu_base - mean(m[,i]) # re-center
}

## ------------------------------------------------------------------------
sharpe <- function(x,rf=0.03/252) {
  sr <- apply(x,2,function(col) {
    er = col - rf
    return(mean(er)/sd(er))
  })
  return(sr)
}

## ------------------------------------------------------------------------
require(pbo)
my_pbo <- pbo(m,s=8,f=sharpe,threshold=0)

## ------------------------------------------------------------------------
summary(my_pbo)

## ------------------------------------------------------------------------
require(lattice)
require(latticeExtra)
require(grid)

histogram(my_pbo,type="density")
xyplot(my_pbo,plotType="degradation")
xyplot(my_pbo,plotType="dominance",increment=0.001)
xyplot(my_pbo,plotType="pairs")
xyplot(my_pbo,plotType="ranks",ylim=c(0,20))
dotplot(my_pbo)

## ----,echo=TRUE,eval=TRUE------------------------------------------------
require(doParallel)

cluster <- makeCluster(2) # or use detectCores()
registerDoParallel(cluster)
p_pbo <- pbo(m,s=8,f=sharpe,allow_parallel=TRUE)
stopCluster(cluster)

summary(p_pbo)

