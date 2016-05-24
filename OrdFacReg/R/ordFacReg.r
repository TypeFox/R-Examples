ordFacReg <- function(D, Z, fact, ordfact, ordering = NA, type = c("LS", "logreg"), intercept = TRUE, display = 0, eps = 0){

if (identical(NA, ordering) == TRUE){ordering <- rep("i", length(ordfact))}

## prepare data
prep <- prepareData(Z, fact, ordfact, ordering, intercept)
Y <- prep$Y
n <- prep$n
p <- prep$p
m <- prep$m
f <- prep$f
c <- p - f
kj <- prep$kj
JJs <- prep$JJs

## generate matrix V for constraints and basis B
mats <- constraintMats(m, c, f, JJs, kj)
B <- mats$B
V <- mats$V
q <- dim(V)[2]

## define functions according to type of regression
## least squares
if (type == "LS"){
    lik.value <- function(beta, D, Y){return(lmSS(beta, D, Y)$L)}
    beta.fct <- function(D, Y, A, JJs, q){return(LSEsubspace(D, Y, A, JJs, q)$beta)}
    dL.fct <- function(beta, D, Y){return(lmSS(beta, D, Y)$dL)}
}

## logistic regression
if (type == "logreg"){
    lik.value <- function(beta, D, Y){return(logRegLoglik(beta, D, Y)$L)}
    beta.fct <- function(D, Y, A, JJs, q){return(logRegSubspace(D, Y, A, JJs, q)$beta)}
    dL.fct <- function(beta, D, Y){return(logRegDeriv(beta, D, Y)$dL)}
}

## start according to Table 2 in Duembgen, Huesler & Rufibach (2007)
A <- 1:q
beta <- beta.fct(D, Y, A, JJs, q)
A <- Abeta(V, beta)
     
## iterative algorithm
iter1 <- 0
crit <- t(B) %*% dL.fct(beta, D, Y)

if (length(A) > 0){
while (max(crit[A]) > eps){

     iter1 <- iter1 + 1
     
     # choose first a such that crit is maximal
     tmp <- (crit[A] == max(crit[A])) * A
     tmp <- tmp[tmp > 0]
     a <- min(tmp)
     A <- A[A != a]     
     beta_new <- beta.fct(D, Y, A, JJs, q)

     iter2 <- 0
     while (max(t(V) %*% beta_new) > eps){ 
          iter2 <- iter2 + 1
          t <- maxStep(beta, beta_new, V, eps)
          beta <- (1 - t) * beta + t * beta_new
          A <- Abeta(V, beta)
          beta_new <- beta.fct(D, Y, A, JJs, q)
          } # end while
     
     beta <- beta_new
     A <- Abeta(V, beta)      
     crit <- t(B) %*% dL.fct(beta, D, Y)
     
     if  (display == 1){print(paste("iter1 = ", iter1, " / L = ", round(lik.value(beta, D, Y), 3), " / A = ", list(A), sep = ""))}          
     if (length(A) == 0){break}
     
} # end while
} 

# re-arrange elements of beta for all factors that are estimated in *decreasing* order
for (u in 1:length(JJs)){if (ordering[u] == "d"){beta[JJs[[u]]] <- beta[rev(JJs[[u]])]}}

## generate output
dimnames(beta) <- list(dimnames(Y)[[2]], NULL)
L <- lik.value(beta, D, Y)
return(list("L" = L, "beta" = beta, "A" = A, "design.matrix" = Y))
}
