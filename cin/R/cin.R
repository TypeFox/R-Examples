cin <- function(X, k=5, type=c("sum", "correlation"), weight=NULL, TR=NULL, interp=FALSE)
{
  type <- match.arg(type, c("sum", "correlation"))
  
  if (is.null(weight)) {
    if (!is.null(TR)) {
      weight <- gammaHRF(TR)
    } else {
      stop("Please specify either TR or weight!")
    }
  } 
  iscor <- !(type=="sum")

  procTest <- function(xatom)  {
    x <- xatom[[1]]
    stc <- xatom[[2]]
    stt <- xatom[[3]]
    
    yc <- truncSum(stc, x, weight, interp, iscor)
    yt <- truncSum(stt, x, weight, interp, iscor)
    intertest(yc, yt, k)
  }

  reMat <- t(sapply(X, procTest))
  ## print(dim(reMat))
  if (nrow(reMat)>1) {
    re <- apply(reMat[,1:3], 2, sum)
    de <- (re[1]-re[2])/sqrt(re[3])
    outlist <- c(re, de, 1-pnorm(de))
  } else {
    outlist <- reMat
  }
  colnames(outlist) <- c("Score","Exp","Var","Dev","p.value")
  outlist <- as.list(outlist)
  class(outlist) <- c("cin")
  return(outlist)
}

