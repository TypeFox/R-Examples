### Initial
library(pbdDEMO, quietly = TRUE)
library(phyclust, quietly = TRUE)

### Load data
data.path <- paste(.libPaths()[1], "/phyclust/data/pony524.phy", sep = "")
pony.524 <- read.phylip(data.path)
X <- pony.524$org
K0 <- 1
Ka <- 2

### Find MLEs
if(comm.rank() == 0){
  ret.K0 <- find.best(X, K0)
  ret.Ka <- find.best(X, Ka)
} else{
  ret.K0 <- NULL
  ret.Ka <- NULL
}

### Broadcast fitted models to all workers
ret.K0 <- bcast(ret.K0)
ret.Ka <- bcast(ret.Ka)
LRT <- -2 * (ret.Ka$logL - ret.K0$logL)

### The user defined function
FUN <- function(jid){
  X.b <- bootstrap.seq.data(ret.K0)$org

  ### find.best is recommended, but more time-consuming 
  # ret.K0 <- find.best(X.b, K0)
  # ret.Ka <- find.best(X.b, Ka)
  ret.K0 <- phyclust(X.b, K0)
  repeat{
    ret.Ka <- phyclust(X.b, Ka)
    if(ret.Ka$logL > ret.K0$logL){
      break
    }
  }

  LRT.b <- -2 * (ret.Ka$logL - ret.K0$logL)
  # cat("comm.rank:", comm.rank(), " jid:", jid, "LRT:", LRT.b, "\n")
  LRT.b
}

### Task pull and summary
ret <- task.pull(1:100, FUN)
if(comm.rank() == 0){
  LRT.B <- unlist(ret) 
  cat("K0: ", K0, "\n",
      "Ka: ", Ka, "\n",
      "logL K0: ", ret.K0$logL, "\n",
      "logL Ka: ", ret.Ka$logL, "\n",
      "LRT: ", LRT, "\n",
      "p-value: ", mean(LRT > LRT.B), "\n", sep = "")
}

### Finish
finalize()
