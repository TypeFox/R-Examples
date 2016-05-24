# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## initializations
library("simFrame")

## start cluster
cl <- makeCluster(2, type="PSOCK")

## load required packages on workers
clusterEvalQ(cl, {
        library("simFrame")
        library("robCompositions")
        library("mvtnorm")
    })

## setup random number stream
clusterSetRNGStream(cl, iseed=12345)

## define function and control class for generating data
crnorm <- function(n, mean, sigma) isomLRinv(rmvnorm(n, mean, sigma))
sigma <- matrix(c(1, -0.5, 1.4, -0.5, 1, -0.6, 1.4, -0.6, 2), 3, 3)
dc <- DataControl(size = 150, distribution = crnorm, 
    dots = list(mean = c(0, 2, 3), sigma = sigma))

## define control class for the insertion of missing values 
nc <- NAControl(NArate = c(0.01, 0.03, 0.05, 0.07, 0.09))

## define function for simulation runs
sim <- function(x, orig) {
    i <- apply(x, 1, function(x) any(is.na(x)))
    ni <- length(which(i))
    xKNNa <- impKNNa(x)$xImp
    xLS <- impCoda(x, method = "lm")$xImp
    c(knn = aDist(xKNNa, orig)/ni, LS = aDist(xLS, orig)/ni)
}

## export object to workers
clusterExport(cl, c("crnorm", "sigma", "dc", "nc", "sim"))

## run simulation
results <- clusterRunSimulation(cl, dc, nrep = 100, NAControl = nc, fun = sim)

## stop cluster
stopCluster(cl)

## inspect results
head(results)
aggregate(results)

## plot results
plot(results, ylab = "Relative Aitchison distance")
simDensityplot(results, NArate=0.07, 
    alpha = 0.6, xlab = "Relative Aitchison distance")
