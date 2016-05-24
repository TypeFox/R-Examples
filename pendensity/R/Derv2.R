Derv2 <- function(penden.env,lambda0) {
  beta.val <- get("beta.val",penden.env)
  b <- length(beta.val)-get("N",penden.env)

  #Derv2.temp is the second order derivative of the likelihood, constructed with the outer product of the first derivative
  Derv2.temp <- matrix(0,b,b)
  for (j in 1:b) {
    Derv2.temp[,j] <- get("Derv1.cal",penden.env)%*%get("Derv1.cal",penden.env)[j,]
  }

  #penalty <- lambda0*get("Dm",penden.env)
  #return(list(Derv2.pen=(Derv2.temp+penalty),Derv2.cal=(Derv2.temp)))
  return(list(Derv2.pen=(Derv2.temp+lambda0*get("Dm",penden.env)),Derv2.cal=(Derv2.temp)))
} 
