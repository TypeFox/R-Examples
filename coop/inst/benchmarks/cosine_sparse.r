library(rbenchmark)
library(slam)
library(memuse)
library(coop)

reps <- 30
cols <- c("test", "replications", "elapsed", "relative")


generate <- function(m, n, size)
{
  x <- matrix(0, m, n)
  x[sample(m*n, size=size)] <- 100
  x
}

m <- 2000
n <- 2000



### Very sparse matrix --- .1% dense
cat("Very sparse (.1% dense):\n")
size <- .001*m*n

dense <- generate(m, n, size)
sparse <- as.simple_triplet_matrix(dense)

memuse(dense)
memuse(sparse)

benchmark(cosine(dense), cosine(sparse), as.matrix(sparse), columns=cols, replications=reps)
cat("\n\n")



### Fairly sparse matrix - 5% dense
cat("Fairly sparse (1% dense):\n")
size <- .01*m*n

dense <- generate(m, n, size)
sparse <- as.simple_triplet_matrix(dense)

memuse(dense)
memuse(sparse)

benchmark(cosine(dense), cosine(sparse), as.matrix(sparse), columns=cols, replications=reps)
