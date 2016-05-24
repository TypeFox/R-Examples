require(inline)
require(Rcpp)

get.model <- function(model.file) {
    do.call(paste,as.list(paste(readLines(model.file),"\n")))
}

src <- get.model("run.model.cpp")

cppbugs.plugin <- getPlugin("RcppArmadillo")
cppbugs.plugin$env$PKG_CXXFLAGS <- "-std=c++0x"
cppbugs.plugin$env$PKG_LIBS <- paste(cppbugs.plugin$env$PKG_LIBS,"-larmadillo")

linear.model <- cxxfunction(signature(XR="numeric", yr="numeric",iterations="integer",burn="integer",adapt="adapt",thin="integer"), includes="#include <cppbugs/cppbugs.hpp>",body=src,settings=cppbugs.plugin)


NR <- 1000
NC <- 5
X <- cbind(rep(1,NR),matrix(rnorm(NR*NC),nrow=NR,ncol=NC))
b <- rnorm(NC + 1)
b[1] <- 10
y <- X %*% b + rnorm(NR)

res <- linear.model(XR=X,yr=y,iterations=1e5L,burn=1e5L,adapt=1e3L,thin=5L)

cat("actual vs estimated:\n")
act.vs.est <- cbind(as.matrix(b),res$b)
colnames(act.vs.est) <- c("actual","estimated")
print(act.vs.est)

cat("ar:",res$ar,"\n")
cat("R^2:",res$rsq,"\n")

