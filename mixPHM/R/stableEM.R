`stableEM` <-
function(x, K, numEMstart = 5, method = "separate", Sdist = "weibull", cutpoint = NULL,
         EMoption = "classification", EMstop = 0.0001, maxiter = 1000, print.likvec = TRUE)
#computes survclust models for different starting solutions and selects the best one
#numEMstart... number of different starting solutions
{

if (is.data.frame(x)) x <- as.matrix(x)

n <- dim(x)[1]
nuvec <- 1:numEMstart
EMst.l <- tapply(nuvec,nuvec, function(y) {                            #list with starting matrices (vectors)
          EMstart <- sample(1:K, n, replace=TRUE)
          return(EMstart)
          })

reslist <- lapply(EMst.l,function(y) {                                #list of models for different EMstart
                 res <- phmclust(x=x,K=K,method=method,Sdist=Sdist,cutpoint = cutpoint, EMstart=y,
                                  EMoption=EMoption,EMstop=EMstop,maxiter=maxiter)
                 return(res)})
likvec <- sapply(reslist,function(y) y$likelihood[length(y$likelihood)])   #likelihoods for all models
if (print.likvec == TRUE) {
     cat("\n")
     cat("Likelihood values for different starting solutions: \n")
     print(likvec)
     }
     
best.ind <- (1:numEMstart)[likvec==max(likvec)]                        #best model
result <- reslist[[best.ind[1]]]
class(result) <- "mws"
result
}

