BoundedAntiMeanTwo <- function(g1, w1, g2, w2, K1 = 1000, K2 = 400, delta = 10 ^ (-4), errorPrec = 10, output = TRUE){

n <- length(g1)

##--------------------------------------------
# generate the G_{s, i}'s
##--------------------------------------------
Gsi <- matrix(NA, nrow = n, ncol = n)
for (s in 1:n){
    c.gw <- cumsum((g1 * w1)[s:n])
    c.w <- cumsum(w1[s:n])
    c.ratio <- rev(c.gw / c.w)
    Gsi[s:n, s] <- rev(cummax(c.ratio))
} # end s
if (output == TRUE){print("Generation of G_{s, i} done")}

##--------------------------------------------
# compute b_n*
##--------------------------------------------
B <- bstar_n(g1, w1, g2, w2)
if (output == TRUE){print("Computation of B done")}

##--------------------------------------------
# initialize the algorithm
##--------------------------------------------
b_old <- rep(B, n)
sg_old <- Subgradient(b_old[-n], g1, w1, g2, w2, B, Gsi)
t_old <- 1 / sqrt(sum((sg_old$grad) ^ 2))
k <- 1
h <- 1
h_add <- 1
error <- Inf

##--------------------------------------------
# run
##--------------------------------------------
while ((k <= K1) & (error > delta)){

    k <- k + 1
    h <- h + h_add
    
    ##--------------------------------------------
    ## compute new candidates
    ##--------------------------------------------
    v_new <- b_old - t_old * sg_old$grad
    b_new <- BoundedAntiMean(y = v_new, w = rep(1 / n, n), a = rep(B, n), b = rep(Inf, n))
    sg_new <- Subgradient(b_new[-n], g1, w1, g2, w2, B, Gsi)
    norm_grad_new <- sqrt(sum((sg_new$grad) ^ 2))
    if(k <= K2){t_new <- 1 / norm_grad_new}
    if(k > K2){t_new <- 1 / (h ^ 0.1 * norm_grad_new)}

    ##--------------------------------------------
    ## re-define quantities
    ##--------------------------------------------
    b_old <- b_new
    t_old <- t_new
    sg_old <- sg_new
    if ((k == 2) | ((k / errorPrec) == round(k / errorPrec))){
        a <- BoundedAntiMean(y = g1, w = w1, a = b_old, b = rep(Inf, n))
        error <- max(abs(b_old - BoundedAntiMean(g2, w2, a = rep(-Inf, n), b = a))) 
        L <- LSfunctional(a, g1, w1, b_old, g2, w2)           
        }
    psi_new <- ""
    if (output == TRUE){print(paste("k = ", k, " / tau = ", disp(t_old, 4), " / Error = ", disp(error, -log10(delta)), " / norm subgradient = ", disp(norm_grad_new, 4), " / LS = ", disp(L, 4), sep = ""))}
}

##--------------------------------------------
# generate output
##--------------------------------------------
a_old <- BoundedAntiMean(y = g1, w = w1, a = b_old, b = rep(Inf, n))
L <- LSfunctional(a_old, g1, w1, b_old, g2, w2)
res <- list("g1" = a_old, "g2" = b_old, "L" = L, "error" = error, "k" = k, "tau" = t_old)
return(res)
}







