library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e1

x <- rstpareto(n, 2)

#test CDF
z <- 0:4/4
ecdf(x)(z)

#test EC
eecf(x)(z)

class(eecf(x))
class(ecdf(x))

print(eecf(x))
print(ecdf(x))

cbind(eecf(x)(sort(x)), 
environment(eecf(x))$"Gx")

print(summary(eecf(x)))
print(summary(ecdf(x)))

plot(eecf(x))
plot(ecdf(x))

