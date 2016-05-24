fun.fit.gl.v3m <-
function(a, b, data, fun, no = 10000, maxmin=TRUE,leap=3,FUN="runif.sobol")
{
e <- a
d <- b - e

if(as.character(substitute(fun)) == "fun.auto.perc.rs") {
param <- "rs"
ncol.init <- 2
init.par <- c(3, 4)
fun1 <- fun.rs.perc.min
g.init <- fun.gen.qrn(n = no, dimension = ncol.init, scrambling = leap,FUN=FUN)* d + e
valid.init <- do.call("cbind", fun.rs.perc.sol.alt(g.init, data = data))
}
if(as.character(substitute(fun)) == "fun.auto.mm.fmkl") {
param <- "fmkl"
ncol.init <- 2
init.par <- c(3, 4)
fun1 <- fun.fmkl.mm.min
g.init <- fun.gen.qrn(n = no, dimension = ncol.init, scrambling = leap,FUN=FUN)* d + e
valid.init <- do.call("cbind", fun.fmkl.mm.sol.alt(g.init, data = data))
}
# Initial values created as above:
# qgl(0,L1,l2,L3,L4,param)<=min(data) && qgl(1,L1,l2,L3,L4,param)>=max(data)
#### work out all the initial values:
### Then: find out which initial values are valid:
test.valid <- fun.check.gld.multi(valid.init[, 1], valid.init[, 2], valid.init[, 3], valid.init[, 4],param = param)
#test.valid <- apply(valid.init, 1, function(x, param)gl.check.lambda.alt(x[1], x[2], x[3], x[4], param = param), param)
valid.init <- valid.init[test.valid,  ]
# This ensures the initial values covers the minimum and the maximum of the data values:
if(maxmin == TRUE) {
maxmin.ind<-fun.minmax.check.gld(data,valid.init,param,0,1)
valid.init <- valid.init[maxmin.ind,  ]
}

# Check to see to make sure that you do not get any Inf or -Inf which will cause the optimisation to fail

optim.check<-is.notinf(optim.fun3.C.m(valid.init,data,param))
valid.init<-valid.init[which(optim.check==TRUE),]

# Work out the values of the minimisation function for these initial values:
test.ind <- apply(valid.init[,3:4], 1, function(x, data, fun1)
fun1(coef = x, data = data), data, fun1)
# Use the best initial values to find the initial solution: The following has been vectorised:
# Work out the minimum:
min.ind <- which(test.ind == (min(test.ind, na.rm = TRUE)))
# Note the following is different, as we now have 4 columns :) init.sol MUST have four columns to allow the values to go through!
init.sol <- matrix(valid.init[min.ind,  ], ncol = 4)
# Now check for existence of the init.sol
if(dim(init.sol)[1] == 0 || is.null(dim(init.sol))) {
stop("Randomized starting points did not result in valid gld")
}
optim.result <- lapply(1:nrow(init.sol), function(i, init.sol, optim.fun3.C, data, param)optim(init.sol[i,  ], optim.fun3.C, data = data, param = param, control = list(maxit = 200000)), init.sol,optim.fun3.C, data, param)
optim.result.m <- matrix(unlist(sapply(1:length(optim.result), function(i, optim.result) optim.result[[i]]$par, optim.result)), ncol = 4, byrow = TRUE)
optim.result.o <- matrix(unlist(sapply(1:length(optim.result), function(i, optim.result)optim.result[[i]]$value, optim.result)), ncol = 1, byrow = TRUE)
# Then do some checks as usual:
# Obtain unique results:
if(dim(optim.result.m)[1] == 1) {
unique.optim.result <- optim.result.m
unique.optim.result.o <- optim.result.o
} else if(dim(optim.result.m)[1] != 1) {
unique.optim.result <- optim.result.m[!duplicated.data.frame(data.frame(signif(optim.result.m,
2))),  ]
unique.optim.result.o <- optim.result.o[!duplicated.data.frame(data.frame(signif(optim.result.o,
2))),  ]
}
# Obtain valid results:
if(is.null(dim(unique.optim.result))) {
unique.optim.result <- matrix(unique.optim.result, nrow = 1)
unique.optim.result.o <- matrix(unique.optim.result.o, nrow = 1)
}
# This step is not necessary at all:
# index<-apply(unique.optim.result,1,function(x,param) gl.check.lambda.alt(x[1],x[2],x[3],x[4],param=param),param)
return(list("unique.optim.result"=unique.optim.result,"unique.optim.result.o"=unique.optim.result.o))
}
