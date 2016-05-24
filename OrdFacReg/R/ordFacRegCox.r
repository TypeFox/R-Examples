ordFacRegCox <- function(ttf, tf, Z, fact, ordfact, ordering = NA, intercept = TRUE, display = 0, eps = 0){

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

## start according to Table 2 in Duembgen, Huesler & Rufibach (2007)
A <- 1:q
beta <- coxSubspace(ttf, tf, Y, A, JJs, q)$beta
A <- Abeta(V, beta)
          
## iterative algorithm
iter1 <- 0
crit <- t(B) %*% coxDeriv(beta, ttf, tf, Y)$dL

if (length(A) > 0){
while ((max(crit[A]) > eps) && (iter1 <= 20)){

     iter1 <- iter1 + 1
     A_old <- A
     
     # choose first a such that crit is minimal
     tmp <- (crit[A] == max(crit[A])) * A
     tmp <- tmp[tmp > 0]
     a <- min(tmp)
     A <- A[A != a]
     beta_new <- coxSubspace(ttf, tf, Y, A, JJs, q)$beta

     iter2 <- 0
     while (max(t(V) %*% beta_new) > eps){ 
          iter2 <- iter2 + 1
          t <- maxStep(beta, beta_new, V, eps)
          beta <- (1 - t) * beta + t * beta_new
          A <- Abeta(V, beta)
          beta_new <- coxSubspace(ttf, tf, Y, A, JJs, q)$beta
          } # end while
     
     beta <- beta_new
     A <- Abeta(V, beta)
     L <- coxLoglik(beta, ttf, tf, Y)$L
     crit <- t(B) %*% coxDeriv(beta, ttf, tf, Y)$dL
     
     if (display == 1){print(paste("iter1 = ", iter1, " / L = ", round(L, 3), " / A = ", list(A), sep = ""))}          
     if (length(A) == 0){break} 
     #if (length(A) == length(A_old)){if (A == A_old){break}}
     
} # end while

} # end if

# re-arrange elements of beta for all factors that are estimated in *decreasing* order
for (u in 1:length(JJs)){if (ordering[u] == "d"){beta[JJs[[u]]] <- beta[rev(JJs[[u]])]}}

## generate output
dimnames(beta) <- list(dimnames(Y)[[2]], NULL)
L <- coxLoglik(beta, ttf, tf, Y)$L
return(list("L" = L, "beta" = beta, "A" = A, "design.matrix" = Y))
}
