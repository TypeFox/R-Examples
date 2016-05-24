"fun.fit.gl.v2a.nw" <-
function(a, b, data, fun, no = 10000, maxmin=TRUE, nclass = nclass.scott,leap=3,FUN="runif.sobol")
{
e <- a
d <- b - e

if(as.character(substitute(fun)) == "fun.auto.perc.rs") {
param <- "rs"
ncol.init <- 2
init.par <- c(3, 4)


fun1 <- fun.rs.perc.min
g.init <- fun.gen.qrn(n = no, dimension = ncol.init, scrambling = leap,FUN=FUN)* d + e

valid.init<-do.call("cbind",fun.rs.perc.sol.alt(g.init, data = data))

#valid.init <- matrix(unlist(apply(g.init, 1, function(x, data)fun.rs.perc.sol(x, data = data), data)), ncol = 4, byrow = TRUE)
}
if(as.character(substitute(fun)) == "fun.auto.mm.fmkl") {
param <- "fmkl"
ncol.init <- 2
init.par <- c(3, 4)
fun1 <- fun.fmkl.mm.min
g.init <- fun.gen.qrn(n = no, dimension = ncol.init, scrambling = leap,FUN=FUN)* d + e
valid.init<-do.call("cbind",fun.fmkl.mm.sol.alt(g.init, data = data))

#valid.init <- matrix(unlist(apply(g.init, 1, function(x, data)
#fun.fmkl.mm.sol(x, data = data), data)), ncol = 4, byrow = TRUE)
}

# Initial values created as above:
# qgl(0,l1,l2,l3,l4,param)<=min(data) && qgl(1,l1,l2,l3,l4,param)>=max(data)
#### work out all the initial values:
### Then: find out which initial values are valid:



test.valid <- gl.check.lambda.alt(valid.init[,1], valid.init[,2], valid.init[,3], valid.init[,4], param = param)

#test.valid <- apply(valid.init, 1, function(x, param)gl.check.lambda.alt(x[1], x[2], x[3], x[4], param = param), param)
valid.init <- valid.init[test.valid,  ]
# This ensures the initial values covers the minimum and the maximum of the data values:
if(maxmin == TRUE) {
test.max <- apply(valid.init, 1, function(x, param)qgl(1, x[1], x[2], x[3], x[4], param = param), param) >= max(data)
test.min <- apply(valid.init, 1, function(x, param)qgl(0, x[1], x[2], x[3], x[4], param = param), param) <= min(data)
valid.init <- valid.init[as.logical(test.max * test.min),  ]
}

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
# Then use this value to optimize a goodness of fit:
if(is.numeric(nclass)) {
breaks <- pretty.su(data, nclass)
} else {
n.class <- nclass(data)
breaks <- pretty.su(data, n.class)
}
data.c <- cut(data, breaks,include.lowest = TRUE)
   counts <- tabulate(data.c, length(levels(data.c)))
   yy<-counts/sum(counts)


optim.result <- lapply(1:nrow(init.sol), function(i, init.sol, optim.fun2.nw, 
data, param, breaks, yy)optim(init.sol[i,  ], optim.fun2.nw, data = data, 
param = param, control = list(maxit = 200000), breaks = breaks, yy = yy), 
init.sol, optim.fun2.nw, data, param, breaks, yy)
optim.result.m <- matrix(unlist(sapply(1:length(optim.result), function(i, 
optim.result)optim.result[[i]]$par, optim.result)), ncol = 4, byrow = TRUE)

optim.result.o <- matrix(unlist(sapply(1:length(optim.result), function(i, 
optim.result)optim.result[[i]]$value, optim.result)), ncol = 1, byrow = TRUE)

# Then do some checks as usual:
# Obtain unique results:
if(dim(optim.result.m)[1] == 1) {
unique.optim.result <- optim.result.m 
unique.optim.result.o <- optim.result.o }

else if(dim(optim.result.m)[1] !=1)

 {
unique.optim.result <- optim.result.m[!duplicated.data.frame(data.frame(signif(optim.result.m,2))),  ]
unique.optim.result.o <- optim.result.o[!duplicated.data.frame(data.frame(signif(optim.result.o,2))),  ]
}

# Obtain valid results:
if(is.null(dim(unique.optim.result))) {
unique.optim.result <- matrix(unique.optim.result, nrow = 1)
unique.optim.result.o <- matrix(unique.optim.result.o, nrow = 1)

}

# This step is not necessary at all:
# index<-apply(unique.optim.result,1,function(x,param) gl.check.lambda.alt(x[1],x[2],x[3],x[4],param=param),param)
return(list("unique.optim.result"=unique.optim.result,"unique.optim.result.o"=unique.optim.result.o))
#return(unique.optim.result)

}

