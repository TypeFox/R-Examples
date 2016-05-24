"fun.fit.gl.v2b.nw" <-
function(a, b, c,d,data, fun, no = 100, maxmin=TRUE, nclass = nclass.scott)
{

if(as.character(substitute(fun)) == "fun.auto.perc.rs") {
param <- "rs"
fun1 <- fun.rs.perc.min
}
if(as.character(substitute(fun)) == "fun.auto.mm.fmkl") {
param <- "fmkl"
fun1 <- fun.fmkl.mm.min
}

init.sol <- c(a,b,c,d)

# Then use this value to optimize a goodness of fit:
if(is.numeric(nclass)) {
breaks <- pretty.su(data, nclass)
} else {
n.class <- nclass(data)
breaks <- pretty.su(data, n.class)
}
data.c <- cut(data, breaks,include.lowest=TRUE)
 counts <- tabulate(data.c, length(levels(data.c)))
   yy<-counts/sum(counts)
optim.result <- 
optim(init.sol, optim.fun2.nw, data = data, param = param, 
control = list(maxit = 200000), breaks = breaks, yy = yy)

return(optim.result$par)
}

