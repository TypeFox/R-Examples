cat("\ntest 1D optimization:")
# Deterministic toy example
fr <- function(v) {   ## Rosenbrock Banana function reduced to 1D
  10 * (1 - v^2)^2 + (1 - v)^2
}
set.seed(123)

# Initial parameter values, including duplicates. See ?init_grid.
parsp <- init_grid(lower=c(x=0),upper=c(x=2))

# add function values
simuls <- cbind(parsp,bb=apply(parsp,1,"fr"))

# optimization
bbresu <- bboptim(simuls)
print(bbresu)

# refine with additional points
while ( ! bbresu$convergence ) {
  candidates <- rbb(bbresu)
  newsimuls <- cbind(candidates,bb=apply(candidates,1,"fr"))
  bbresu <- bboptim(rbind(bbresu$fit$data,newsimuls))
}

expect_equal(bbresu$optr$par,1,tolerance=1e-4)
