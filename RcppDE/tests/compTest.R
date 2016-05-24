
Wild <- function(x) { 		## 'Wild' function, global minimum at about -15.81515
    sum(10 * sin(0.3 * x) * sin(1.3 * x^2) + 0.00001 * x^4 + 0.2 * x + 80)/length(x)
}

Rastrigin <- function(x) {
    sum(x+2 - 10 * cos(2*pi*x)) + 20
}

Genrose <- function(x) { 	## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}


maxIt <- 25                        # not excessive but so that we get some run-time on simple problems

suppressMessages(library(DEoptim)) 	# the original, currently 2.0.7
suppressMessages(library(RcppDE))   # the contender

basicDE <- function(n, maxIt, fun) DEoptim::DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n),
                                                    control=list(NP=10*n, itermax=maxIt, trace=FALSE))#, bs=TRUE))
cppDE <- function(n, maxIt, fun) RcppDE::DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n),
                                                 control=list(NP=10*n, itermax=maxIt, trace=FALSE))#, bs=TRUE))

set.seed(42)
valBasic <- basicDE(5, maxIt, function(...) Rastrigin(...))
set.seed(42)
valCpp <- cppDE(5, maxIt, function(...) Rastrigin(...))
#stopifnot( all.equal(valBasic, valCpp) )

runPair <- function(n, maxIt, fun, funname) {
    gc()
    set.seed(42)
    bt <- system.time(invisible(ores <- basicDE(n, maxIt, fun)))[3]

    gc()
    set.seed(42)
    xptr <- .Call("putFunPtrInXPtr", funname, package="RcppDE")
    ct <- system.time(invisible(cres <- cppDE(n, maxIt, xptr)))[3]

    #stopifnot(all.equal(ores, cres))

    gc()
    set.seed(42)
    rt <- system.time(invisible(rres <- cppDE(n, maxIt, fun)))[3]

    #stopifnot(all.equal(ores, rres))

    return(data.frame(DEoptim=bt, RcppDEc=ct, RcppDEr=rt))
}

cat("# At", format(Sys.time()), "\n")

reps <- c(5, 10, 20, 50)

res <- rbind(do.call(rbind, lapply(reps, runPair, maxIt, function(...) Rastrigin(...), "rastrigin")),
             do.call(rbind, lapply(reps, runPair, maxIt, function(...) Wild(...), "wild")),
             do.call(rbind, lapply(reps, runPair, maxIt, function(...) Genrose(...), "genrose"))
             )
res <- rbind(res, colMeans(res))

rownames(res) <- c(paste("Rastrigin", reps, sep=""),
                   paste("Wild", reps, sep=""),
                   paste("Genrose", reps, sep=""),
                   "MEANS")

res$ratioRcppCompToBasic <- res[,2]/res[,1]
res$pctGainOfRcppComp <- (1-res[,2]/res[,1])*100
#res$netSpeedUpC <- res[,1]/res[,2]

res$ratioRcppRToBasic <- res[,3]/res[,1]
res$pctGainOfRcppR <- round((1-res[,3]/res[,1])*100, digits=3)
#res$netSpeedUpR <- res[,1]/res[,3]


print(res)
cat("# Done", format(Sys.time()), "\n")
