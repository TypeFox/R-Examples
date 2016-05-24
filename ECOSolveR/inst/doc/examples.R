## ----echo=F--------------------------------------------------------------
### get knitr just the way we like it

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)


## ------------------------------------------------------------------------
library(ECOSolveR)
library(Matrix)
set.seed(182391)
n <- 1000L
m <- 10L
density <- 0.01
c <- c(rep(0.0, n), rep(1.0, n))

## ------------------------------------------------------------------------
sprandn <- function(nrow, ncol, density) {
    items <- ceiling(nrow * ncol * density)
    matrix(c(rnorm(items),
             rep(0, nrow * ncol - items)),
           nrow = nrow)
}

## ------------------------------------------------------------------------
A <- sprandn(m, n, density)
Atilde <- Matrix(cbind(A, matrix(rep(0.0, m * n), nrow = m)), sparse = TRUE)
b <- rnorm(m)
I <- diag(n)
G <- rbind(cbind(I, -I),
           cbind(-I, -I))
G <- Matrix(G, sparse = TRUE)
h <- rep(0.0, 2L * n)
dims <- list(l = 2L * n, q = NULL, e = 0L)

## ------------------------------------------------------------------------
## Solve the problem
z <- ECOS_csolve(c = c, G = G, h = h, dims = dims, A = Atilde, b = b)

## ------------------------------------------------------------------------
names(z)
z$infostring

## ------------------------------------------------------------------------
x <- z$x[1:n]
u <- z$x[(n+1):(2*n)]
nnzx = sum(abs(x) > 1e-8)
sprintf("x reconstructed with %d non-zero entries", nnzx / length(x) * 100)

