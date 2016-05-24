bayestwoNorm <- function(xx,yy,R=10^4){

# now fit the distributions to the data using bayesm


X <- matrix(1,length(xx),1)
Bbar <- c(0,0)
A <- as.matrix(10^-3)


b <- matrix(double(R*2),ncol=2)
S <- b

Data1 <- list(y=as.matrix(xx,ncol=1),x=X)
Data2 <- list(y=as.matrix(yy,ncol=1),x=X)

Prior <- list(betabar=Bbar,A)
Mcmc <- list(R=R)


out1 <- runireg(Data=list(y=xx,X=X), Prior=list(betabar=Bbar[1],A=A), Mcmc=list(R=R))
out2 <- runireg(Data=list(y=yy,X=X), Prior=list(betabar=Bbar[2],A=A), Mcmc=list(R=R))

b[,1] <- out1$betadraw
b[,2] <- out2$betadraw

S[,1] <- out1$sigmasqdraw
S[,2] <- out2$sigmasqdraw



model <- list()
model$b <- b
model$S <- S
return(model)


}