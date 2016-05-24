library(pbdDEMO, quietly=TRUE)

init.grid()

comm.set.seed(1234, diff=TRUE)

# default 2x2 grid
dx.0 <- ddmatrix("rnorm", nrow=500, ncol=200, bldim=2)
print(dx.0)

# 1x4 grid
dx.1 <- redistribute(dx.0, ICTXT=1)
print(dx.1)

# 4x1 grid
dx.2 <- redistribute(dx.1, ICTXT=2)
print(dx.2)

# default 2x2 grid with different blocking
dx_new <- redistribute(dx.2, bldim=c(200,200), ICTXT=0)
print(dx_new)

finalize()
