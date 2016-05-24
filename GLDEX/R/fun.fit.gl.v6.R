`fun.fit.gl.v6` <-
function (a, b, data, fun, no = 10000, maxmin = TRUE, 
    leap = 3, FUN = "runif.sobol",trial.n=100,len=1000,type=7) 
{
    e <- a
    d <- b - e
    if (as.character(substitute(fun)) == "fun.auto.perc.rs") {
        param <- "rs"
        ncol.init <- 2
        init.par <- c(3, 4)
        fun1 <- fun.rs.perc.min
        g.init <- fun.gen.qrn(n = no, dimension = ncol.init, 
            scrambling = leap, FUN = FUN) * d + e
        valid.init <- do.call("cbind", fun.rs.perc.sol.alt(g.init, 
            data = data))
    }
    if (as.character(substitute(fun)) == "fun.auto.mm.fmkl") {
        param <- "fmkl"
        ncol.init <- 2
        init.par <- c(3, 4)
        fun1 <- fun.fmkl.mm.min
        g.init <- fun.gen.qrn(n = no, dimension = ncol.init, 
            scrambling = leap, FUN = FUN) * d + e
        valid.init <- do.call("cbind", fun.fmkl.mm.sol.alt(g.init, 
            data = data))
    }
    test.valid <- gl.check.lambda.alt(valid.init[, 1], valid.init[, 
        2], valid.init[, 3], valid.init[, 4], param = param)
    valid.init <- valid.init[test.valid,]
    if (maxmin == TRUE) {
        test.max <- apply(valid.init, 1, function(x, param) qgl(1, 
            x[1], x[2], x[3], x[4], param = param), param) >= 
            max(data)
        test.min <- apply(valid.init, 1, function(x, param) qgl(0, 
            x[1], x[2], x[3], x[4], param = param), param) <= 
            min(data)
        valid.init <- valid.init[as.logical(test.max * test.min), 
            ]
    }
    optim.check <- is.notinf(apply(valid.init, 1, function(x, 
        data, param,trial.n,type) optim.fun6(x, data, param,trial.n,type), data, 
        param,trial.n,type))
    valid.init <- valid.init[optim.check,,drop=F]
    if(is.matrix(valid.init)==FALSE){
    return("Failed to obtain a suitable set of initial values")
    }
    test.ind <- apply(valid.init[, 3:4,drop=F], 1, function(x, data, 
        fun1) fun1(coef = x, data = data), data, fun1)
    min.ind <- which(test.ind == (min(test.ind, na.rm = TRUE)))
    init.sol <- matrix(valid.init[min.ind, ], ncol = 4)
    if (dim(init.sol)[1] == 0 || is.null(dim(init.sol))) {
        stop("Randomized starting points did not result in valid gld")
    }
    optim.result <- lapply(1:nrow(init.sol), function(i, init.sol, 
        optim.fun6, data, param,len,type) optim(init.sol[i, ], optim.fun6, 
        data = data, param = param, len=len,type=type,
        control = list(maxit = 2e+05)), init.sol, optim.fun6, data, 
        param,len,type)
    optim.result.m <- matrix(unlist(sapply(1:length(optim.result), 
        function(i, optim.result) optim.result[[i]]$par, optim.result)), 
        ncol = 4, byrow = TRUE)
    optim.result.o <- matrix(unlist(sapply(1:length(optim.result), 
        function(i, optim.result) optim.result[[i]]$value, optim.result)), 
        ncol = 1, byrow = TRUE)
    if (dim(optim.result.m)[1] == 1) {
        unique.optim.result <- optim.result.m
        unique.optim.result.o <- optim.result.o
    }
    else if (dim(optim.result.m)[1] != 1) {
        unique.optim.result <- optim.result.m[!duplicated.data.frame(data.frame(signif(optim.result.m, 
            2))), ]
        unique.optim.result.o <- optim.result.o[!duplicated.data.frame(data.frame(signif(optim.result.o, 
            2))), ]
    }
    if (is.null(dim(unique.optim.result))) {
        unique.optim.result <- matrix(unique.optim.result, nrow = 1)
        unique.optim.result.o <- matrix(unique.optim.result.o, 
            nrow = 1)
    }
return(list("unique.optim.result"=unique.optim.result,
"unique.optim.result.o"=unique.optim.result.o))
}

