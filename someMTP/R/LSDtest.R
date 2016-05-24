lsd.test <- function (resp, alternative = 1, null = NULL, D = NULL, data=NULL){

  call <- match.call()
  data=.getXY(resp,alternative,null,data); rm(resp,alternative,null)
  n <- dim(data$Y)[1]
  p <- dim(data$Y)[2]
  k <- dim(data$X)[2]
  
  if((!is.null(data$Z)&&(dim(data$Z)[2]>0)) ) {
    h <- dim(data$Z)[2]
    IP0 <- diag(n) - data$Z%*%solve(t(data$Z)%*%data$Z)%*%t(data$Z)
	} else{
    h <- 0
	IP0 <- diag(n) 
   }

  if (is.null(D)) {
   D=diag(t(data$Y)%*%IP0%*%data$Y)  
  } else   if(is.function(D)) D <- D(resp=data$Y,null=data$Z)
  D <- as.matrix(D)  
  q <- dim(D)[2]
   
  H = t(data$Y)%*%IP0%*%data$X%*%solve(t(data$X)%*%IP0%*%data$X)%*%t(data$X)%*%IP0%*%data$Y
  G <- t(data$Y)%*%IP0%*%data$Y - H
  DHD	 <- t(D) %*% H %*% D
  DGD <- t(D) %*% G %*% D
  
  
  out <- new("lsd.object")  
  out @ call = call
  out @ df = c(q-1+k,n-h-k+1-q)
  out @ F = out@df[2]/out@df[1]*sum(diag(DHD))/sum(diag(DGD))
  out @ globalP = 1-pf(out@F,out@df[1],out@df[2])
  out @ D = D
  rownames(out @ D) = colnames(data$Y)
  
  return(out)
}