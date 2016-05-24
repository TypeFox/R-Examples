## ------------------------------------------------------------------------
library(ldamatch)
set.seed(257)

SIZE <- 15
condition <- as.factor(c(rep("control", 2 * SIZE), rep("treatment", SIZE)))
covariate1 <- c(rnorm(2 * SIZE), rnorm(SIZE, 1, 2))

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariate1, t_halt)
print(table(condition, is.in))

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariate1, t_halt, method = "monte")
print(table(condition, is.in))

## ------------------------------------------------------------------------
covariate2 <- c(rnorm(2 * SIZE), rnorm(SIZE, 1, 2))
covariates <- cbind(covariate1, covariate2)

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariates, t_halt)
print(table(condition, is.in))

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariates, t_halt, method = "monte")
print(table(condition, is.in))

## ------------------------------------------------------------------------
my.props <- prop.table(c(control = 4, treatment = 3))
is.in <- ldamatch(condition, covariates, U_halt, props = my.props)
print(table(condition, is.in))

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariates, wilks_halt)
print(table(condition, is.in))

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariates, wilks_halt, method = "monte")
print(table(condition, is.in))

## ------------------------------------------------------------------------
is.in <- ldamatch(condition, covariates, ad_halt, method = "heuristic")
print(table(condition, is.in))

## ------------------------------------------------------------------------
combined_halting_test <- create_halting_test(c(t_halt, ad_halt))
threshes <- c(.2, .02)
is.in <- ldamatch(condition, covariates, combined_halting_test, threshes)
print(table(condition, is.in))

## ------------------------------------------------------------------------
estimate_exhaustive(min_preserved = 42, condition, cases_per_second = 100)
foreach::registerDoSEQ()
is.ins <- ldamatch(condition, covariate1, t_halt, method = "exhaustive")
print(table(condition, is.ins[[1]]))
print(length(is.ins))

# (Confirm exhaustive search by applying heuristic search to it.)
is.in <- ldamatch(condition[is.ins[[1]]], covariate1[is.ins[[1]]], t_halt)
print(table(condition[is.ins[[1]]], is.in))

