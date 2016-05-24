DearBeggMonotone <- function(y, u, lam = 2, maxiter = 1000, CR = 0.9, NP = NA, trace = TRUE){

# lam: first weight (either 1 or 2, see Dear & Begg, p. 239. They recommend lam = 2.)

# general parameters
n <- length(y)
k <- 1 + floor(n / 2)
teststat <- abs(y) / u
p <- 2 * pnorm(-teststat)
p0 <- p

## data preparation: sort all vectors in decreasing order of p-values
ind <- order(p) 
ind <- rev(ind)
p <- p[ind]
y <- y[ind]
u <- u[ind]
teststat <- teststat[ind]

if (is.na(NP)){size <- 10 * (k + 2)} else {size <- NP}
inipop <- matrix(runif(size * (k + 2)), ncol = k + 2, nrow = size, byrow = TRUE)
for (i in 1:nrow(inipop)){inipop[i, ] <- c(sort(inipop[i, 1:k]), runif(1, -20, 20), runif(1, 0, 20))}

d0 <- DEoptim::DEoptim(fn = DearBeggToMinimize, lower = c(rep(0, k), -20, 0), upper = c(rep(1, k), 20, 50), 
    control = DEoptim.control(strategy = 2, bs = FALSE, NP = size, trace = trace, itermax = maxiter, CR = CR, F = 0.8, 
    initialpop = inipop), y, u, lam) 

w <- as.numeric(d0$optim$bestmem)[1:k] 
theta <- as.numeric(d0$optim$bestmem)[k + 1] 
sigma <- as.numeric(d0$optim$bestmem)[k + 2] 
hij <- Hij(theta, sigma, y, u, teststat)$Hij  
ll.DE <- DearBeggLoglik(w, theta, sigma, y, u, hij, lam)$LL

res <- list("w" = w, "theta" = theta, "sigma" = sigma, "p" = p, "y" = y, "u" = u, "loglik" = ll.DE, "DEoptim.res" = d0)
return(res)
}
