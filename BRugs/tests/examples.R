source("local.R")

if (test) {

library("BRugs")

## Selected examples which take a few seconds in total to run

test.models <- c("Air", "Asia", "Beetles", "BiRats", "Camel",
                 "Dugongs", "Dyes", "Equiv", "Eyes",
                 "Line", "OtreesMVN", "Rats", "Stacks",
                 "Surgical", "Surgicalrand")

test.params <- list(Air = c("X", "theta"),
                    Asia = c("bronchitis", "either", "lung.cancer"),
                    Beetles = c("alpha", "beta", "rhat"),
                    BiRats = c("mu.beta", "sigma"),
                    Camel = c("Sigma2", "rho", "tau"),
                    Dugongs = c("U3","alpha", "beta", "gamma", "sigma"),
                    Dyes = c("sigma2.btw", "sigma2.with", "theta"),
                    Equiv = c("equiv", "mu", "phi", "pi","sigma1", "sigma2", "theta"),
                    Eyes = c("P", "lambda", "sigma"),
                    Line = c("alpha", "beta", "sigma"),
                    OtreesMVN = c("mu","sigma", "sigmaC"),
                    Rats = c("alpha0", "beta.c", "sigma"),
                    Stacks = c("b", "b0", "outlier[21]","outlier[3]", "outlier[4]", "sigma"),
                    Surgical = "p",
                    Surgicalrand = c("p","pop.mean", "sigma")
                    )

test.modelfile <- paste(test.models,"model.txt",sep="")
test.datafile <- paste(test.models,"data.txt",sep="")
test.inits <- paste(test.models,"inits.txt",sep="")
test.pattern <- paste("^", test.models, ".*\\.txt$", sep="")
### Test for posterior means within 10 percent of previously saved values

res.true <- dget(file="examples.stats.R")

exfiles <- unlist(lapply(test.pattern, function(tp) dir(options()$OpenBUGSExamples, pattern=tp, full.names=TRUE)))
ok <- file.copy(unique(exfiles), tempdir())
if(!all(ok)) 
    stop("Some files could not be copied from OpenBUGS examples to the temporary directory")

for (i in seq(along=test.models)) {
    fit <- BRugsFit(data=test.datafile[i], inits=test.inits[i], 
                    modelFile=test.modelfile[i], para=test.params[[test.models[i]]],
                    nBurnin=5000, nIter=20000, nThin=1, numChains=1, seed=1,
                    working.directory=tempdir())
    stopifnot(isTRUE(all.equal(fit$Stats$mean, res.true[[i]]$Stats$mean, tol=1e-01)))
}

}
