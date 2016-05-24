context("Check bivariate GLM case")
require(MASS)
#Examples of "sim" 
set.seed (1)
J <- 15
n <- J*(J+1)/2
group <- rep (1:J, 1:J)
mu.a <- 5
sigma.a <- 2
a <- rnorm (J, mu.a, sigma.a)
b <- -3
x <- rnorm (n, 2, 1)
sigma.y <- 6
y <- rnorm (n, a[group] + b*x, sigma.y)
u <- runif (J, 0, 3)
y123.dat <- cbind (y, x, group)
# Linear regression 
x1 <- y123.dat[,2]
y1 <- y123.dat[,1]
M1 <- glm (y1 ~ x1)

cases <- data.frame(x1 = seq(-2, 2, by=0.1))
sim.results <- gelmansim(mod = M1, newdata = cases, n.sims=200, na.omit=TRUE)


test_that("returned dataframe is correct size", {
  expect_that(dim(sim.results)[1], equals(length(seq(-2, 2, by=0.1))))
  expect_that(dim(sim.results)[2], equals(4))
  expect_that(sim.results, is_a("data.frame"))
})

test_that("values of simulations are sensible", {  
  expect_that(all(sim.results$yhatMin < sim.results$yhats), is_true())
  expect_that(all(sim.results$yhatMax > sim.results$yhats), is_true())
  expect_that(all(!is.na(sim.results$yhats)), is_true())
  expect_that(all(!is.na(sim.results$yhatMin)), is_true())
  expect_that(all(!is.na(sim.results$yhatMax)), is_true())
})

context("Check multivariate GLM case")

dat <- as.data.frame(y123.dat)
M2 <- glm (y1 ~ x1 + group, data=dat)

cases <- expand.grid(x1 = seq(-2, 2, by=0.1), 
                    group=seq(1, 14, by=2))

sim.results <- gelmansim(M2, newdata=cases, n.sims=200, na.omit=TRUE)

test_that("returned dataframe is correct size", {
  expect_that(dim(sim.results)[1], equals(nrow(cases)))
  expect_that(dim(sim.results)[2], equals(3 + ncol(cases)))
  expect_that(sim.results, is_a("data.frame"))
})

test_that("values of simulations are sensible", {  
  expect_that(all(sim.results$yhatMin < sim.results$yhats), is_true())
  expect_that(all(sim.results$yhatMax > sim.results$yhats), is_true())
  expect_that(all(!is.na(sim.results$yhats)), is_true())
  expect_that(all(!is.na(sim.results$yhatMin)), is_true())
  expect_that(all(!is.na(sim.results$yhatMax)), is_true())
})

context("Check multivariate GLM case with factor")

dat$group <- factor(dat$group)
M3 <- glm (y1 ~ x1 + group, data=dat)

cases <- expand.grid(x1 = seq(-2, 2, by=0.1), 
                     group=seq(1, 14, by=2))
cases$group <- factor(cases$group)

sim.results <- gelmansim(M3, newdata=cases, n.sims=200, na.omit=TRUE)

test_that("returned dataframe is correct size", {
  expect_that(dim(sim.results)[1], equals(nrow(cases)))
  expect_that(dim(sim.results)[2], equals(3 + ncol(cases)))
  expect_that(sim.results, is_a("data.frame"))
})

test_that("values of simulations are sensible", {  
  expect_that(all(sim.results$yhatMin < sim.results$yhats), is_true())
  expect_that(all(sim.results$yhatMax > sim.results$yhats), is_true())
  expect_that(all(!is.na(sim.results$yhats)), is_true())
  expect_that(all(!is.na(sim.results$yhatMin)), is_true())
  expect_that(all(!is.na(sim.results$yhatMax)), is_true())
})

context("Check multivariate LM case with factor")

dat$group <- factor(dat$group)
M4 <- lm(y1 ~ x1 + group, data=dat)

cases <- expand.grid(x1 = seq(-2, 2, by=0.1), 
                     group=seq(1, 14, by=2))
cases$group <- factor(cases$group)

sim.results <- gelmansim(M4, newdata=cases, n.sims=200, na.omit=TRUE)


test_that("returned dataframe is correct size", {
  expect_that(dim(sim.results)[1], equals(nrow(cases)))
  expect_that(dim(sim.results)[2], equals(3 + ncol(cases)))
  expect_that(sim.results, is_a("data.frame"))
})

test_that("values of simulations are sensible", {  
  expect_that(all(sim.results$yhatMin < sim.results$yhats), is_true())
  expect_that(all(sim.results$yhatMax > sim.results$yhats), is_true())
  expect_that(all(!is.na(sim.results$yhats)), is_true())
  expect_that(all(!is.na(sim.results$yhatMin)), is_true())
  expect_that(all(!is.na(sim.results$yhatMax)), is_true())
})
