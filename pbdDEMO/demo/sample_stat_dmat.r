### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial MPI.
library(pbdDEMO, quietly = TRUE)
init.grid()
if(comm.size() != 2){
  comm.stop("This example requries 2 processors.")
}

comm.set.seed(1234, diff = TRUE)

# --------------------------
# set up distributed matrix
# --------------------------

n <- 25
p <- 15

mean <- 100
sd <- 25

bldim <- c(2,2)

#dx <- ddmatrix("rnorm", nrow=n, ncol=p, bldim=bldim, mean=mean, sd=sd)
x <- matrix(rnorm(n*p), n, p)
dx <- as.ddmatrix(x)

# --------------------------
# sample statistics
# --------------------------

# means
mn <- mean(dx)
comm.print(mn)

rm <- rowMeans(dx)
comm.print(rm)

cm <- colMeans(dx)
comm.print(cm)

# variances
sd1 <- apply(dx, MARGIN=2, FUN=stats::sd, reduce=T)
sd2 <- sd(dx, reduce=T) # unlike R, no complaint for doing this
                         # in fact, it's much faster than the above
comm.print(sd1)
comm.print(sd2)

cv_diag <- diag( cov(dx)[1:5, 1:5] )
comm.print(cv_diag)


finalize()
