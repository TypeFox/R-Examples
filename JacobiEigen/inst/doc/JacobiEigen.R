## ---- comment = ""-------------------------------------------------------
suppressMessages(library(dplyr))
library(JacobiEigen)
library(stats)

imod <- aov(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, iris)
(R <- cor(resid(imod)))

rEig <- JacobiR(R) 
cEig <- Jacobi(R)
identical(rEig, cEig)  ## the R and Rcpp implementations are identical
cEig
(eEig <- eigen(R))
all.equal(eEig$values, cEig$values)  ## eigenvalues are (practically) identical
crossprod(eEig$vectors, cEig$vectors) %>% ## eigenvectors differ in signs
  round(10) 

## ---- comment = ""-------------------------------------------------------
library(microbenchmark)
microbenchmark(JacobiR(R), Jacobi(R), eigen(R))

## ---- comment = ""-------------------------------------------------------
suppressMessages(library(tidyr))
set.seed(1234)
N <- 100
iseq <- seq(5, 50, by = 5)
res <- lapply(iseq,  function(n) {
  S <- crossprod(matrix(rnorm(N*n), N, n))/N
  runTime <- microbenchmark(JacobiR(S), Jacobi(S), eigen(S), times = 20)
  c(n = n, with(runTime, tapply(time, expr, median))/1000)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame %>% 
  gather(key = expr, value = time, `JacobiR(S)`, `Jacobi(S)`, `eigen(S)`)

suppressMessages(library(ggplot2))
ggplot(res) + aes(x = n, y = log10(time), colour = expr) + geom_line() + geom_point() +
  theme(legend.position = "top", legend.title = element_blank()) + xlab("matrix size") +
  ylab(expression(log[10]("median run time in milliseconds")))

## ---- echo=FALSE, results="asis"-----------------------------------------
cat("\\newpage\n")

## ---- echo=FALSE, results="asis"-----------------------------------------
cat("\\newpage\n")

