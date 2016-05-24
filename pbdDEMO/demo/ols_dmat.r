### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial MPI.
library(pbdDEMO, quietly = TRUE)
init.grid()
if(comm.size() != 4){
  comm.stop("This example requries 4 processors.")
}


n <- 1250
p <- 40

mean <- 100
sd <- 1000

ymin <- 0
ymax <- 500

bldim <- c(2,2)

comm.set.seed(1234, diff=TRUE)

dx <- ddmatrix("rnorm", nrow=n, ncol=p, bldim=bldim, mean=mean, sd=sd)
dy <- ddmatrix("runif", nrow=n, ncol=1, bldim=bldim, min=ymin, max=ymax)

mdl <- lm.fit(dx, dy)

dx.new <- ddmatrix("rnorm", nrow=1, ncol=p, bldim=bldim, mean=mean, sd=sd)
pred <- dx.new %*% mdl$coefficients

comm.cat(paste("\nThe predicted y value is:", submatrix(pred), "\n"), quiet=T)

finalize()
