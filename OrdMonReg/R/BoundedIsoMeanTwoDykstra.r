BoundedIsoMeanTwoDykstra <- function(g1, w1, g2, w2, K1 = 1000, delta = 10 ^ (-8), output = TRUE){

n <- length(g1)
g <- c(g2, g1)
w <- c(w2, w1)

r <- 3   # number of cones
error <- 1
Ii <- matrix(0, ncol = 2 * n, nrow = r)    # note that I_i is a vector of same dimension as g!
iter <- 1

while ((iter <= K1) & (error > delta)){

    # n modulo r according to definition on p. 838 in Dykstra (1983)
    nmodr <- iter %% r
    if (nmodr == 0){nmodr <- r}
    if (nmodr == 1){cone <- minK1}
    if (nmodr == 2){cone <- minK2}
    if (nmodr == 3){cone <- minK3}

    g.new <- cone(g - Ii[nmodr, ], w, n)
    Ii[nmodr, ] <- g.new - (g - Ii[nmodr, ])

    error <- sum((g.new - g) ^ 2)
    g <- g.new

    iter <- iter + 1
    if ((iter > r) & (identical(output, TRUE))){
        L <- LSfunctional(g[(n + 1):(2 * n)], g1, w1, g[1:n], g2, w2)
        print(paste("iteration = ", iter, " / error = ", disp(error, -log10(delta)), " / LS = ", disp(L, 4), sep = ""))}
}

## generate output
a <- g[(n + 1):(2 * n)]
b <- g[1:n]
L <- LSfunctional(a, g1, w1, b, g2, w2)
    
res <- list("g1" = a, "g2" = b, "L" = L, "error" = error, "k" = iter)
return(res)
}
