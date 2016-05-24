HMT <- function(N = 5000, H = 10){
   x <- rgamma(n = N, shape = 2, scale = 5)
   p <- 0.04 * (x^(-3/2)) * ((8 + 5*x)^2)
   lam <- (8 + 5*x)/(1.25 * (x^(3/2)))
   y <- rgamma(n = N, shape = p, rate = lam)
   pop <- cbind(x,y)
   pop <- pop[order(pop[, 1]), ]

   brk <- (0:H) * sum(pop[, 1]/H)
   brk[length(brk)] <- cumsum(pop[, 1])[N]
   strat <- cut(cumsum(pop[, 1]), brk)

   pop <- cbind(strat, pop)
   pop
}
