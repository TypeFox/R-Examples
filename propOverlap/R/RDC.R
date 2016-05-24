RDC <- function(GMask, Y){
  Y <- as.factor(Y)
  levels(Y) = c(1,2)
  Relative.DC <- max.col(cbind(rowSums(GMask[,Y == 1])/sum(Y == 1), rowSums(GMask[,Y == 2])/sum(Y == 2)))
  return(Relative.DC)
}