## ------------------------------------------------------------------------
set.seed(2)
J <- 500
mu <- - 2*abs(rnorm(J))
frac_sig <- 0.1
X <- rbinom(J, 1, frac_sig)
t1 <- rnorm(J, X*mu,1)
t2 <- rnorm(J, X*mu,1)

## ---- fig.show='hold'----------------------------------------------------
scatter.smooth(t1,t2)

## ---- fig.show='hold'----------------------------------------------------
P_current <- pnorm(t2)

## ------------------------------------------------------------------------
alpha <- 0.05
P_adjusted <- p.adjust(P_current,"bonferroni")
ind <- (P_adjusted<alpha)
cat(c("Number of significant tests using Bonferroni: ", sum(ind) ))


## ------------------------------------------------------------------------
q <- alpha/J
sigma <- rep(1,J)
source("../R/bayes_weights.R")
res <- bayes_weights(t1,sigma,q)
w <- res$w
P_weighted <- P_current/w
P_w_adjusted <- p.adjust(P_weighted,"bonferroni")
ind_w <- (P_w_adjusted<alpha)
cat(c("Number of significant tests using Weighting: ", sum(ind_w) ))

## ------------------------------------------------------------------------
which(ind==1)
which(ind_w==1)

## ------------------------------------------------------------------------
x = 266
P_current[x]
t2[x]

## ------------------------------------------------------------------------
w[x]
t1[x]

## ------------------------------------------------------------------------
plot(t1,w)

## ------------------------------------------------------------------------
P_prior <- pnorm(t1)

## ------------------------------------------------------------------------
N_current <- 1
N_prior <- 1

## ------------------------------------------------------------------------
source("../R/iGWAS.R")
res_unw <- iGWAS(P_current, N_current, P_prior, N_prior, weighting_method="unweighted")

## ------------------------------------------------------------------------
res_w <- iGWAS(P_current, N_current, P_prior, N_prior, sides=1)

