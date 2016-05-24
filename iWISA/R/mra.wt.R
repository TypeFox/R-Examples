mra.wt <-
function(x.wt)
{
wf<-attr(x.wt,"wavelet")
J<-length(x.wt)-1
method<-attr(x.wt,"class")
boundary<-attr(x.wt,"boundary")
if(method=="modwt") n<-length(x.wt[[1]])
else n<-2*length(x.wt[[1]])

x.mra <- vector("list", J + 1)
zero <- vector("list", J + 1)
names(zero) <- c(paste("d", 1:J, sep = ""), paste("s", J,
sep = ""))
class(zero) <- method
attr(zero, "wavelet") <- wf
attr(zero, "boundary") <- boundary
zero[[J + 1]] <- x.wt[[J + 1]]
if (method == "modwt") {
for (k in 1:J) zero[[k]] <- numeric(n)
x.mra[[J + 1]] <- imodwt(zero)
}
else {
for (k in 1:J) zero[[k]] <- numeric(n/2^k)
x.mra[[J + 1]] <- idwt(zero)
}
for (j in J:1) {
zero <- vector("list", j + 1)
names(zero) <- c(paste("d", 1:j, sep = ""), paste("s",
j, sep = ""))
class(zero) <- method
attr(zero, "wavelet") <- wf
attr(zero, "boundary") <- boundary
zero[[j]] <- x.wt[[j]]
if (method == "modwt") {
if (j != 1) {
for (k in c(j + 1, (j - 1):1)) zero[[k]] <- numeric(n)
}
else {
zero[[j + 1]] <- numeric(n)
}
x.mra[[j]] <- imodwt(zero)
}
else {
zero[[j + 1]] <- numeric(n/2^j)
if (j != 1) {
for (k in (j - 1):1) zero[[k]] <- numeric(n/2^k)
}
x.mra[[j]] <- idwt(zero)
}
}
names(x.mra) <- c(paste("D", 1:J, sep = ""), paste("S", J,
sep = ""))
if (boundary == "reflection") {
for (j in (J + 1):1) x.mra[[j]] <- x.mra[[j]][1:(n/2)]
return(x.mra)
}
else {
return(x.mra)
}
}
