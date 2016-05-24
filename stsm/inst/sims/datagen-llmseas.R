#
# Generate a ensemble from the local level plus seasonal model
# The data set 'llmseas' available in the 'stsm' package
# is generated as shown below
#

library("stsm")

# parameters of the data generating process

dgp.n <- 120
dgp.var1 <- 300
dgp.var2 <- 10
dgp.var3 <- 100

dgp.s <- 4
dgp.n0 <- 20

# number of series to generate and seed

iter <- 1000
seed <- 123

# generate data

m <- stsm.model(model = "llm+seas", y = ts(seq_len(dgp.n), frequency = dgp.s), 
  pars = c(var1 = dgp.var1, var2 = dgp.var2, var3 = dgp.var3))

ss <- char2numeric(m, FALSE)
SigmaEV <- eigen(ss$Q)

My <- matrix(nrow = dgp.n, ncol = iter)

set.seed(seed)

for (i in seq_len(iter))
{
  My[,i] <- datagen.stsm(n = dgp.n, 
    model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q), 
    SigmaEV = SigmaEV, n0 = dgp.n0, freq = dgp.s,
    old.version = TRUE)$data
}

llmseas <- ts(My, frequency = dgp.s)
