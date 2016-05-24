library(HyperbolicDist)

### Create vector of nu values and of x values
nus <- c(0:5, 10, 20)
x <- seq(1, 4, length.out = 11)
### Specify the difference between the numerator and denominator orders
k <- 3

### Matrix for results of calculations
results <- matrix(nrow = length(nus), ncol = length(x))
for (i in 1:length(nus)){
    for (j in 1:length(x)) {
        raw <- besselK(x[j], nus[i] + k)/besselK(x[j], nus[i])
        scaled <- besselK(x[j], nus[i] + k, expon.scaled = TRUE)/
                  besselK(x[j], nus[i], expon.scaled = TRUE)
        results[i,j] <- raw/scaled
    }
}

### Should yield 1 in every case
results
max(abs(results - 1))

### Try a negative value for the difference of the orders
k <- -3

results <- matrix(nrow = length(nus), ncol = length(x))
for (i in 1:length(nus)){
    for (j in 1:length(x)) {
        raw <- besselK(x[j], nus[i] + k)/besselK(x[j], nus[i])
        scaled <- besselK(x[j], nus[i] + k, expon.scaled = TRUE)/
                  besselK(x[j], nus[i], expon.scaled = TRUE)
        results[i,j] <- raw/scaled
    }
}

results
max(abs(results - 1))

### Now use besselRatio function
nus <- c(0:5, 10, 20)
x <- seq(1, 4, length.out = 11)
k <- 3

raw <- matrix(nrow = length(nus), ncol = length(x))
scaled <- matrix(nrow = length(nus), ncol = length(x))
compare <- matrix(nrow = length(nus), ncol = length(x))

for (i in 1:length(nus)){
    for (j in 1:length(x)) {
        raw[i,j] <- besselRatio(x[j], nus[i],
                                orderDiff = k)
        scaled[i,j] <- besselRatio(x[j], nus[i],
                                orderDiff = k, useExpScaled = 1)
        compare[i,j] <- raw[i,j]/scaled[i,j]
    }
}
raw
scaled
compare
