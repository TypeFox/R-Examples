
# parameters of the data generating process

dgp.n <- 120
dgp.var1 <- 1600
dgp.var2 <- 100

dgp.n0 <- 20

# number of series to generate and seed

iter <- 1000
seed <- 123

# generate data

set.seed(seed)

My <- matrix(nrow = dgp.n, ncol = iter)
for (i in seq(iter))
  My[,i] <- cumsum(rnorm(dgp.n, sd = sqrt(dgp.var2))) + 
    rnorm(dgp.n, sd = sqrt(dgp.var1))

llm <- ts(My, frequency = 1)
