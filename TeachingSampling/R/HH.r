HH <- function(y,pk){
  y <- as.data.frame(y)
  m <- length(pk)
  Total <- matrix(NA,nrow=3,ncol=dim(y)[2])
  rownames(Total)=c("Estimation", "Standard Error","CVE")
  colnames(Total) <- names(y)
  
  for(k in 1:dim(y)[2]){
    ty <- sum(y[,k]/pk)/m
    Vty <- (1/m)*(1/(m-1))*sum((y[,k]/pk-ty)^2)
    CVe <- 100*sqrt(Vty)/ty 
    Total[,k] <- c(ty,sqrt(Vty),CVe)
  }
  return(Total)
}