recodeData <- function(X){
    NCOL <- ncol(X)
    NROW <- nrow(X)
    result <- matrix(NA,ncol=NCOL,nrow=NROW)
    rownames(result) <- rownames(X)
    colnames(result) <- colnames(X)
    X[X=="0"] <- "XX"
    for(i in 1:NCOL){
      tempVal <- X[,i]
      tempGeno <- columnLabels(tempVal)
      result[,i] <- as.numeric(factor(tempVal,tempGeno))
    }
    result - 1
} 
