prepareData <- function(Z, fact = NA, ordfact, ordering = NA, intercept = TRUE){

if (identical(NA, ordering) == TRUE){ordering <- rep("i", length(ordfact))}

Z <- as.matrix(Z)
if (intercept == TRUE){
    Z <- cbind("intercept" = 1, Z)
    fact <- fact + 1
    ordfact <- ordfact + 1
    }

ordfact.orig <- ordfact
    
if (is.na(max(fact)) == 0){

    # compute number of levels per unordered factor
    u.f <- length(fact)     #unordered factors
    u.kj <- NULL
    for (l in 1:u.f){u.kj[l] <- length(apply(as.matrix(Z[, fact[l]]), 2, table))}

    # add dummys for unordered factors
    Y <- cbind(Z[, -c(fact, ordfact)])
    for (l in 1:u.f){Y <- cbind(Y, dummy(Z[, fact[l]], dimnames(Z)[[2]][fact[l]]))}
    
    # recompute ordfact
    ordfact <- length(Y[1, ]) + (1:length(ordfact))
    Z <- cbind(Y, Z[, ordfact.orig])
}

n <- dim(Z)[1]
p <- dim(Z)[2]
f <- length(ordfact)   #ordered factors

# compute number of levels per ordered factor
kj <- NULL
for (l in 1:f){kj[l] <- length(apply(as.matrix(Z[, ordfact[l]]), 2, table))}

# compute column indices in matrix Y, per factor
JJs <- indexDummy(kj, f, p)

# generate Y, the matrix containing the dummy variables
Y <- cbind(Z[, -ordfact])
for (l in 1:f){Y <- cbind(Y, dummy(Z[, ordfact[l]], dimnames(Z)[[2]][ordfact[l]]))}

# re-arrange columns of Y for all factors that are sought to be estimated in *decreasing* order
for (u in 1:length(JJs)){if (ordering[u] == "d"){Y[, JJs[[u]]] <- Y[, rev(JJs[[u]])]}}

# new dimension
m <- dim(Y)[2]

return(list("Z" = Z, "Y" = Y, "ordfact" = ordfact, "n" = n, "p" = p, "m" = m, "f" = f, "kj" = kj, "JJs" = JJs))
}
