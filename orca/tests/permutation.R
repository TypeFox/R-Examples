library(orca)
data(usastates)
n <- max(usastates)
p <- sample(1:n)
usastates.p <- apply(usastates, c(1,2), function(x) p[x])

o <- count5(usastates)
op <- count5(usastates.p)
all(o == op[p,])

e <- ecount5(usastates)
ep <- ecount5(usastates.p)
all(e == ep)
