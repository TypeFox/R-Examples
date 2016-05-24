best.unbal.haar <-
function(x, criterion = inner.prod.max) {
n <- length(x)
if (n < 2) stop("Input vector too short")
tree <- list(matrix(0, 5, 1))
tree[[1]][1,1] <- 1
ipi <- inner.prod.iter(x)
ind.max <- criterion(x)
tree[[1]][3,1] <- 1
tree[[1]][4,1] <- ind.max
tree[[1]][5,1] <- n
tree[[1]][2,1] <- ipi[ind.max]

j <- 1

while(sum(tree[[j]][5,] - tree[[j]][3,] - rep(1, dim(tree[[j]])[2]))) {
tree <- c(tree, list(matrix(0, 5, 0)))
no.parent.coeffs <- dim(tree[[j]])[2]
no.child.coeffs <- 0
for (i in 1:no.parent.coeffs) {
if (tree[[j]][4,i] - tree[[j]][3,i] >= 1) {
no.child.coeffs <- no.child.coeffs + 1
tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, no.child.coeffs)
tree[[j+1]][1,no.child.coeffs] <- 2*tree[[j]][1,i]-1
ipi <- inner.prod.iter(x[(tree[[j]][3,i]):(tree[[j]][4,i])])
ind.max <- criterion(x[(tree[[j]][3,i]):(tree[[j]][4,i])])
tree[[j+1]][2,no.child.coeffs] <- ipi[ind.max]
tree[[j+1]][3,no.child.coeffs] <- tree[[j]][3,i]
tree[[j+1]][5,no.child.coeffs] <- tree[[j]][4,i]
tree[[j+1]][4,no.child.coeffs] <- ind.max + tree[[j]][3,i] - 1
}
if (tree[[j]][5,i] - tree[[j]][4,i] >= 2) {
no.child.coeffs <- no.child.coeffs + 1
tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, no.child.coeffs)
tree[[j+1]][1,no.child.coeffs] <- 2*tree[[j]][1,i]
ipi <- inner.prod.iter(x[(tree[[j]][4,i] + 1):(tree[[j]][5,i])])
ind.max <- criterion(x[(tree[[j]][4,i] + 1):(tree[[j]][5,i])])
tree[[j+1]][2,no.child.coeffs] <- ipi[ind.max]
tree[[j+1]][3,no.child.coeffs] <- tree[[j]][4,i] + 1
tree[[j+1]][5,no.child.coeffs] <- tree[[j]][5,i]
tree[[j+1]][4,no.child.coeffs] <- ind.max + tree[[j]][4,i]
}
}
j <- j + 1
}

smooth <- sum(x) / sqrt(n)

z <- list(tree=tree, smooth=smooth)


return(z)
}

