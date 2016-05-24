cutMatrix <- function(A, mode="col")
{

if(!is.matrix(A))
  stop("'A' is not of type matrix.")

if(!is.character(mode)){
  cat("'mode' has to be of type character. Default is used. \n")
  mode <- "col"
}
if (!(mode %in% c("col", "row"))){
  cat("'mode='", mode, "'' is not supported. Default is used. \n", sep="")
  mode <- "col"
}

A <- A-min(A)
dimOrg <- dim(A)

if(mode=="col"){
  Ay <- ncol(A)
  n1 <- 1
  n2 <- Ay
  icolS <- which(colSums(A)<=1e-10)
  if(length(icolS)!=0){
    L1 <- length(icolS)
    index1 <- 1:L1
    index2 <- Ay:(Ay-L1+1)
    icolS2 <- icolS[L1:1]
    icolMin <- which(index1==icolS)
    if (length(icolMin)!=0)
      icolMin <- max(icolMin)+1
    else
      icolMin <- 1
    icolMax <- which(index2==icolS2)
    if (length(icolMax)!=0)
      icolMax <- Ay - max(icolMax)
    else
      icolMax <- Ay
    A <- A[,icolMin:icolMax]
    pat <- c(icolMin,icolMax)
  }

} else if(mode=="row"){
  Ax <- nrow(A)
  n1 <- 1
  n2 <- Ax
  irowS <- which(rowSums(A)<=1e-10)
  if(length(irowS)!=0){
    L1 <- length(irowS)
    index1 <- 1:L1
    index2 <- Ax:(Ax-L1+1)
    irowS2 <- irowS[L1:1]
    irowMin <- which(index1==irowS)
    if (length(irowMin)!=0)
      irowMin <- max(irowMin)+1
    else
      irowMin <- 1
    irowMax <- which(index2==irowS2)
    if (length(irowMax)!=0)
      irowMax <- Ax - max(irowMax)
    else
      irowMax <- Ax
    A <- A[irowMin:irowMax,]
    pat <- c(irowMin,irowMax)
  }

}

return(list(A=A, dimOrg=dimOrg, pattern=pat))

}