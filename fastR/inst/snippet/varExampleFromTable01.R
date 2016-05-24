y <- 1:4
prob <- c(0.05,0.20,0.40,0.35)
mean.y <- sum(y*prob); mean.y              # E(Y)
sum((y-mean.y)^2 * prob)                   # Var(Y)
sum(y^2 *prob) - mean.y^2                  # Var(Y) again
