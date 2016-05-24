
Delta <- function(tree, i, lambda1, lambda2){
m <- nrow(tree$merge)
if (m > i){
j <- m  - ancestor(tree, i) + 1
i <- m  - i + 1
spec <- smaller.clade.spectrum2(tree)
l0 <- logratio(spec[j,1],spec[j,2],lambda1, lambda2)
l1 <- logratio(spec[i,1],spec[i,2],lambda1, lambda2)
return(lst =list(delta = l0 - l1, clade.s = spec[c(j,i),]) )
}
else
stop("i must be a non-root ancestral node label")
}

