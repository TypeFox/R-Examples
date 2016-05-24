# a = shape; b = scale (acts like 1/lambda from exponential dist)

a <- 2; b <- 3
m <- b * gamma(1 + 1/a); m                            # mean
v <- b^2 * ( gamma(1 + 2/a) - gamma(1+1/a)^2 ); v     # var
qweibull(0.5,a,b)                                     # median
pweibull(m,a,b)                                       # less than mean
pweibull(6,a,b) - pweibull(1.5,a,b)                   # in a range
pweibull(m+sqrt(v),a,b) - pweibull(m - sqrt(v) ,a,b)  # near mean

# with roles of parameters reversed 

a <- 3; b <- 2
m <- b * gamma(1 + 1/a); m                            # mean
v <- b^2 * ( gamma(1 + 2/a) - gamma(1+1/a)^2 ); v     # var
qweibull(0.5,a,b)                                     # median
pweibull(m,a,b)                                       # less than mean
pweibull(6,a,b) - pweibull(1.5,a,b)                   # in a range
pweibull(m+sqrt(v),a,b) - pweibull(m - sqrt(v) ,a,b)  # near mean

