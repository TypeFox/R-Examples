"errormatrix" <-
function(true, predicted, relative=FALSE)
# if `relative=TRUE', rows are nomalized by their sum,
# the last row is is divided by the total number of misclassifications,
# and the lower right cell is the total misclassification rate.
{
  stopifnot(length(true)==length(predicted))
  tnames <- 
    if(is.factor(true)) levels(true) 
    else unique(true)
  pnames <- 
    if(is.factor(predicted)) levels(predicted) 
    else unique(predicted)
  allnames <- sort(union(tnames, pnames))
  n <- length(allnames)
  true <- factor(true, levels = allnames)
  predicted <- factor(predicted, levels = allnames)
  tab <- table(true, predicted)
  mt <- tab * (matrix(1, ncol = n, nrow = n) - diag( , n, n))
  rowsum <- rowSums(mt)
  colsum <- colSums(mt)
  result <- rbind(cbind(tab, rowsum), c(colsum, sum(colsum)))
  dimnames(result) <- list("true" = c(allnames, "-SUM-"), 
    "predicted" = c(allnames, "-SUM-"))
  if(relative){
    total <- sum(result[1:n, 1:n])
    # normalize last row:
    n1 <- n + 1
    result[n1, 1:n] <- 
        if(result[n1, n1] != 0) result[n1, 1:n] / result[n1, n1] 
        else 0
    # normalize remaining matrix:
    rownorm <- function(Row,Length)
    { return( if(any(Row[1:Length]>0)) Row/sum(Row[1:Length])
              else rep(0,Length+1) ) 
    }
    result[1:n,] <- t(apply(result[1:n,], 1, rownorm, Length=n))
    # normalize lower right cell:
    result[n1, n1] <- result[n1, n1] / total
  }
  return(result)
}
