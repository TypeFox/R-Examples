library(pcalg)

## Tests whether "intervention matrices" can correctly be mapped to a 
## targets/target index pair

nreps <- 30
p <- 20
n <- 200

for (i in 1:nreps) {
  set.seed(i)
  A <- matrix(as.logical(rbinom(n*p, 1, 0.01)), nrow = n, ncol = p)
  
  target.list <- mat2targets(A)
  if (any(duplicated(target.list$targets)))
    stop("Targets are not unique!")
  for (j in 1:n)
    if (!all.equal(which(A[j, ]), 
        target.list$targets[[target.list$target.index[j]]]))
      stop("Targets not correctly represented!")
}

