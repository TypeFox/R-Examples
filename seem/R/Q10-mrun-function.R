# Q10 sensitivity

Q10.mruns <- function(Temp, param.nom, Topt.sens, Tmax.sens, q10.sens) {

nT <- length(Temp)
mat<- matrix(1:4,2,2,byrow=T)
layout(mat,c(3.5,3.5),c(3.5,3.5),respect=TRUE)
par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

# loop to change Topt, Tamx and q10
for(j in 1:3){

# determine parameter
 if(j==1) p <- Topt.sens
 if(j==2) p <- Tmax.sens
 if(j==3) p <- q10.sens

# dimension arrays according to how many temp and param values
 np <- length(p); ncalc <- nT*np
 X  <- matrix(nrow=nT,ncol=np) 

# reset to nominal
 Topt<- param.nom[1]; Tmax <- param.nom[2]; q10 <- param.nom[3]

# set varying parameter and calculate
for (i in 1:np){
  if(j==1) Topt <- p[i]
  if(j==2) Tmax <- p[i]
  if(j==3) q10 <- p[i]
  X[,i] <- Q10(Temp, param=c(Topt, Tmax, q10))[,2]  
}

leg.par <- c("Topt", "Tmax","q10")
nom.par <- c(rep(paste("q10=",param.nom[3]),2),"")
matplot(Temp, X, type="l", lty=1:np, col=1, ylab="Temp effect", xlab="Temp(deg C)")
legend(0,0.5, lty=1:np,  col=1, legend=paste(leg.par[j],"=",p))
mtext(side=1,line=-2,text=nom.par[j])
 if(j==1) X.Topt <- X
 if(j==2) X.Tmax <- X
 if(j==3) X.q10 <- X

} # end of j loop

return(list(X.Topt=X.Topt, X.Tmax=X.Tmax, X.q10=X.q10))

} # end function
