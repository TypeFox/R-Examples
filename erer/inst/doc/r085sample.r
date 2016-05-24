# Sample from a distribution
set.seed(10); rnorm(n = 3, mean = 0, sd = 1)
set.seed(20); runif(n = 4, min = 0, max = 1)

# Random sampling on a vector
set.seed(40); ma <- sample(x = 11:15, size = 3)
set.seed(40); mb <- sample(x = 11:15, size = 8, replace = TRUE); ma; mb

# Random sampling on a data frame
bob <- data.frame(inc = 1:5, year = 2001:2005)
set.seed(40); sam <- sample(x = 1:nrow(bob), size = nrow(bob) - 2)
bob2 <- bob[sam, ]
bob; sam; bob2