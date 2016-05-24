bayesMVN <- function(x,y,R=10^4){

# now fit the distributions to the data using bayesm

Y <- cbind(x,y)
X <- matrix(1,length(x),1)
Bbar <- c(0,0)
A <- 10^-3
nu <- 2
V <- 2*diag(2)

b <- matrix(double(R*2),ncol=2)
S <- matrix(double(R*4),ncol=4)

for (i in 1:R){
  out <- rmultireg(Y,X,Bbar,A,nu,V)
  b[i,] <- out$B
  S[i,] <- out$Sigma
}


model <- list()
model$b <- b
model$S <- S   
return(model)


}