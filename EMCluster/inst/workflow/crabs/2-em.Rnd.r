rm(list = ls(all.names = TRUE))

library(EMCluster)
library(MASS)

x <- as.matrix(crabs[, 4:8])
p <- ncol(x)
min.n <- p * (p + 1) / 2
min.n.iter <- 10
.EMC$short.iter <- 5000

ret <- list()
ret.save <- list()
k <- 1
i.rep <- 0
llhdval.curr <- -Inf
K.max <- 7
repeat{
  # Initial.
  i.rep <- i.rep + 1
  seed <- 1234 + k + i.rep
  set.seed(seed)

  ret[[k]] <- list()
  ret[[k]]$seed <- seed

  # Run RndEM.
  method <- "Rnd.EM"
  repeat{
    tmp.init <- init.EM(x, nclass = k, min.n = min.n, min.n.iter = min.n.iter,
                        method = method)
    #if(any(tmp.init$nc <= min.n)){
    #  next
    #}
    var <- LTSigma2variance(tmp.init$LTSigma)
    flag <- 0
    for(k.var in 1:dim(var)[3]){
      # Check degenerate.
      tmp <- try(solve(var[,, k.var]), silent = TRUE)
      if(class(tmp) == "try-error"){
        flag <- 1
        break
      }
    }
    if(flag == 0){
      ret[[k]]$Rnd <- tmp.init
      break
    }
  }

  # Check llhdval.
  flag <- 1
  if(ret[[k]]$Rnd$llhdval > llhdval.curr){
    llhdval.curr <- ret[[k]]$Rnd$llhdval
    ret.save[[k]] <- ret[[k]]$Rnd
    ret.save[[k]]$case.save <- "Rnd"
    flag <- 0
  } 
  if(flag == 0){
    cat("k = ", k, ", llhdval = ", llhdval.curr, ", ", sep = "")
    cat("nc =", ret.save[[k]]$nc, "\n")
    k <- k + 1
  }

  if(k > K.max){
    break
  }
}
save(list = c("x", "ret", "ret.save"), file = "./data/crabs.RndEM.rda")
