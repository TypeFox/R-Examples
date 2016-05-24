SPA <-
function(Y,# Phenotypes
 G,# Genotypes
 g.size, # Size of the two genes in the gene pair
 cutoff=0.05, # Define a p-value cutoff for the overall effect test
 par=NULL,# Initial values for variance components under null model of interaction test
 est.alt=FALSE # logical; should variance components be estimated under alternative model
)
{
n <- length(Y)
L1 <- g.size[1]; L2 <- g.size[2]
Gene1 <- G[, 1:L1]; Gene2 <- G[, (L1+1):ncol(G)]
w1 <- rep(1, L1)
w2 <- rep(1, L2)
K1 <- KERNEL(Gene1, w1)
K2 <- KERNEL(Gene2, w2)
K3 <- K1*K2
test_o <- Score.Test.Overall(Y, K1, K2, K3)
#test_int <- est_res <- NULL
if (test_o$p.value<cutoff)
{
if (is.null(par))
{
grid <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
test_i <- est <- vector("list", length(grid))
for (i in 1 : length(grid))
{
initials <- c(var(Y), rep(grid[i], 2))
test_i[[i]] <- Score.Test.Interact(Y, K1, K2, K3, initials, method="BFGS", test=TRUE)
}
}

if (!is.null(par))
{
initials <- par
test_i <- list(Score.Test.Interact(Y, K1, K2, K3, initials, method="BFGS", test=TRUE))
}
test.lr <- c()
for (i in 1 : length(test_i))
{
test.lr[i] <- test_i[[i]]$restricted.logLik
}
test_int <-test_i[[which.max(test.lr)]]
if (est.alt)
{
initials <- c(test_int$VCs, 0)
est_res <- Score.Test.Interact(Y, K1, K2, K3, initials, method="BFGS", test=FALSE)
res <- list(test.overall=test_o, 
      test.interaction=test_int, 
    parameter.est.alter=est_res)
} else
{res <- list(test.overall=test_o, 
       test.interaction=test_int)}

} else
{res <- list(test.overall=test_o)}

return(res)
}
