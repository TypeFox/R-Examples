source("helper-diversitree.R")

context("Early Burst")

data(bird.orders)
set.seed(1)
x <- structure(rnorm(length(bird.orders$tip.label)),
               names=bird.orders$tip.label)

p0 <- c(.1, 0)

## With the VCV approach
## Note here that the value of theta has no effect on the
## likelihood.  This is a bit surprising, actually.
lik.vcv <- make.eb(bird.orders, x, control=list(method="vcv"))
fit1 <- find.mle(lik.vcv, p0)

lik.pruning <- make.eb(bird.orders, x, control=list(method="pruning"))
fit2 <- find.mle(lik.pruning, p0)

# This is generating an error in a way that bm and ou don't, but it's
# in R's functions, so might be an OSX/valgrind error again.
ign <- make.bm(bird.orders, x,
               control=list(method="pruning", backend="C"))
ign <- make.ou(bird.orders, x,
               control=list(method="pruning", backend="C"))
lik.pruning.C <- make.eb(bird.orders, x,
                         control=list(method="pruning", backend="C"))
fit3 <- find.mle(lik.pruning.C, p0)

## Fits all behave somewhat sensibly:
expect_that(fit1$lnLik, equals(fit3$lnLik, tolerance=1e-6))
expect_that(fit2$lnLik, equals(fit3$lnLik, tolerance=1e-6))
expect_that(coef(fit1), equals(coef(fit2), tolerance=0.002))
expect_that(coef(fit2), equals(coef(fit3), tolerance=0.002))

## Compare against geiger:
set.seed(1)
fit4 <- suppressWarnings(fitContinuous(bird.orders, x, model="EB",
                                       control=list(niter=5)))
expect_that(fit4$opt$lnL, equals(fit1$lnLik, tolerance=1e-6))
p.g <- coef(fit4)[2:1]
names(p.g)[1] <- "s2"
expect_that(p.g, equals(coef(fit1), tolerance=0.002))

expect_that(lik.vcv(p.g), equals(fit4$opt$lnL))
expect_that(lik.pruning(p.g), equals(fit4$opt$lnL, tolerance=1e-6))
expect_that(lik.pruning.C(p.g), equals(fit4$opt$lnL))

lik.geiger <- fit4$lik

p2 <- c(0.02, -0.1)
l2 <- as.numeric(lik.geiger(p2[2:1]))

expect_that(lik.vcv(p2), equals(l2))
expect_that(lik.pruning(p2), equals(l2))
expect_that(lik.pruning.C(p2), equals(l2))

## These are quite different -- more than I'd expect.
expect_that(lik.pruning(p2), equals(lik.pruning.C(p2)))
