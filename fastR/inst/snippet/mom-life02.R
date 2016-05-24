m <- mean(life01); m
n <- length(life01); n
v <- var(life01) * (n-1) / n; v
#alpha
m^2/v
#lambda
m/v
