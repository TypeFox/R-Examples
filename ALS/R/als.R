"als" <- function(CList, PsiList, S=matrix(), WList=list(), thresh =
                  .001, maxiter = 100, forcemaxiter=FALSE, optS1st=TRUE,
                  x=1:nrow(CList[[1]]), x2 = 1:nrow(S), baseline=FALSE,
                  fixed=vector("list", length(PsiList)), uniC=FALSE,
                  uniS=FALSE, nonnegC = TRUE, nonnegS = TRUE,
                  normS=0, closureC=list())
{
  RD <- 10^20
  PsiAll <- do.call("rbind", PsiList)
  resid <- vector("list", length(PsiList))
  # if a weighting specification is absent, set weights to unity
  if(length(WList) == 0){
    WList <- vector("list", length(PsiList))
    for(i in 1:length(PsiList))
      WList[[i]] <- matrix(1, nrow(PsiList[[1]]), ncol(PsiList[[1]]))
  }
  W <- do.call("rbind",WList)
  # initialize residual matrices to zero
  for(i in 1:length(PsiList)) resid[[i]] <- matrix(0, nrow(PsiList[[i]]),
                                                   ncol(PsiList[[i]]))
  # determine the residuals at the starting values
  for(j in 1:length(PsiList)) {
    for(i in 1:nrow(PsiList[[j]])) {
      resid[[j]][i,] <- PsiList[[j]][i,] - CList[[j]][i,] %*% t(S * WList[[j]][i,])
    }
  }
  # set the initial residual sum of squares 
  initialrss <- oldrss <- sum(unlist(resid)^2)
  cat("Initial RSS", initialrss, "\n")
  iter <- 1
  b <- if(optS1st) 1 else 0
  oneMore <- FALSE
  
  while( ((RD > thresh || forcemaxiter ) && maxiter >= iter) || oneMore) {
    if(iter %% 2 == b)   ## solve for S, get RSS 
      S <- getS(CList, PsiAll, S, W, baseline, uniS, nonnegS, normS, x2)
    else   ## solve for CList, get resid matrices in this step only 
      CList <- getCList(S, PsiList, CList, WList, resid, x, baseline,
                         fixed, uniC, nonnegC, closureC)
    # determine the residuals 
    for(j in 1:length(PsiList)) {
      for(i in 1:nrow(PsiList[[j]])) {
        resid[[j]][i,] <- PsiList[[j]][i,] - CList[[j]][i,] %*% t(S * WList[[j]][i,])
      }
    }
    rss <- sum(unlist(resid)^2)
    RD <- ((oldrss - rss) / oldrss)
    oldrss <- rss
    typ <- if(iter %% 2 == b) "S" else "C"
    cat("Iteration (opt. ", typ, "): ", iter, ", RSS: ", rss, ", RD: ", RD,
        "\n", sep = "")
    iter <- iter + 1
    ## make sure the last iteration enforces any normalization/closure
    oneMore <- (normS > S && (iter %% 2 != b) && maxiter != 1) ||
    (length(closureC) > 0 && (iter %% 2 == b) ) 
  }
  cat("Initial RSS / Final RSS =", initialrss, "/", rss, "=",
      initialrss/rss,"\n")
  return(list(CList = CList, S = S, rss = rss, resid = resid, iter = iter))
}

