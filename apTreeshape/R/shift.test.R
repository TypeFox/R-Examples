shift.test <- function(tree, node, lambda1 = 1, lambda2 = 100, nrep = 1000, silent = FALSE){
if (class(tree)[1] != "treeshape") stop("Error: object 'tree' not of class 'treeshape'")
if ((node+2) > length(tree$names)) stop("Error: node > number of taxa - 2 ")
obj <- Delta(tree, node, lambda1, lambda2)
if (obj$clade.s[1,1]%%2 == 0) 
{L <- obj$clade.s[1,1]/2  
l0 <- sapply(sample(1:L, nrep, replace = TRUE),FUN = function(l){logratio(obj$clade.s[1,1],l,lambda1, lambda2)})
}
else 
{
L <- (obj$clade.s[1,1]+1)/2
l0 <- sapply(sample(1:L, nrep, replace = TRUE,prob = c(rep(2,L-1),1)),FUN = function(l){logratio(obj$clade.s[1,1],l,lambda1, lambda2)})
}
if (obj$clade.s[2,1] > 3){
if (obj$clade.s[2,1]%%2 == 0) 
{L <- obj$clade.s[2,1]/2  
l1 <- sapply(sample(1:L, nrep, replace = TRUE),FUN = function(l){logratio(obj$clade.s[2,1],l,lambda1, lambda2)})
}
else 
{
L <- (obj$clade.s[2,1]+1)/2
l1 <- sapply(sample(1:L, nrep, replace = TRUE,prob = c(rep(2,L-1),1)),
FUN = function(l){logratio(obj$clade.s[2,1],l,lambda1, lambda2)})
}}
else
{
l1 <- logratio(obj$clade.s[2,1],1,lambda1, lambda2)
l1 <- rep(l1, nrep)
}
Pleft <- mean( (l0-l1) <= obj$delta ) 
Pright <- mean((l0-l1) >= obj$delta ) 
Delta <- obj$delta
P <- min(Pleft, Pright)
if (silent) return(P) 
else{
cat("Test of diversification rate shift at node", node, "\n")
cat("Delta1 statistic = ", Delta, "\n")
cat("P-value = ", P, "\n")
cat("")
cat("alternative hypothesis: the diversification rate shift is less than ", lambda2/lambda1," fold 
the ancestral rate \n")
cat("Note: The P-value was computed using ", nrep," Monte-Carlo replicates. \n")}
}
