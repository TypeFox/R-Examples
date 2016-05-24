
demo.LargeBenchmark  <- function() {

    Wild <- function(x) { 		## 'Wild' function, global minimum at about -15.81515
        sum(10 * sin(0.3 * x) * sin(1.3 * x^2) + 0.00001 * x^4 + 0.2 * x + 80)/length(x)
    }

    Rastrigin <- function(x) {
        sum(x+2 - 10 * cos(2*pi*x)) + 20
    }

    Genrose <- function(x) { 		## One generalization of the Rosenbrock banana valley function (n parameters)
        n <- length(x)
        1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
    }

    maxIt <- 250                        # not excessive but so that we get some run-time on simple problems

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
    stopifnot( all.equal(valBasic, valCpp) )

    runPair <- function(n, maxIt, fun) {
        gc()
        set.seed(42)
        bt <- system.time(invisible(basicDE(n, maxIt, fun)))[3]

        gc()
        set.seed(42)
        ct <- system.time(invisible(cppDE(n, maxIt, fun)))[3]

        return(data.frame(DEoptim=bt, RcppDE=ct))
    }

    cat("# At", format(Sys.time()), "\n")

    reps <- c(50, 100, 200)

    res <- rbind(do.call(rbind, lapply(reps, runPair, maxIt, function(...) Rastrigin(...))),
                 do.call(rbind, lapply(reps, runPair, maxIt, function(...) Wild(...))),
                 do.call(rbind, lapply(reps, runPair, maxIt, function(...) Genrose(...)))
                 )
    res <- rbind(res, colMeans(res))

    rownames(res) <- c(paste("Rastrigin", reps, sep=""),
                       paste("Wild", reps, sep=""),
                       paste("Genrose", reps, sep=""),
                       "MEANS")

    res$ratioRcppToBasic <- res[,2]/res[,1]
    res$pctGainOfRcpp <- (1-res[,2]/res[,1])*100
    res$netSpeedUp <- res[,1]/res[,2]

    print(res)
    cat("# Done", format(Sys.time()), "\n")
}

demo.LargeBenchmark()
