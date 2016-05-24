# exampleOF.R -- version 2010-12-28
nR <- 100L; nC <- 5L 
X <- array(rnorm(nR * nC),dim = c(nR,nC)); y <- rnorm(nC)

OF1 <- function(param, X, y) max(abs(y - X %*% param))

createF <- function(X, y) {
    newFun <- function(param) max(abs(y - X %*% param))
    newFun
}
OF2 <- createF(X, y)

param <- rnorm(nC)
OF1(param, X, y)
OF2(param)

remove(list = c("X","y"), envir = .GlobalEnv); ls()
OF1(param,X,y)
OF2(param)

environment(OF2)
whereToLook <- environment(OF2)
ls(envir = whereToLook)