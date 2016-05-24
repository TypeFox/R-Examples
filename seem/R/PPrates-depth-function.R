# PPrates-depth-function.R

# smith depth averaged
PPT.Smith <- function(Ls,k,z,Pmax,alpha){
 y <- (Pmax/alpha)^2
 L <- Ls*exp(-k*z)
 if(Ls > 0) PPTSmith <- (Pmax/(k*z)) * log( (Ls+sqrt(y+Ls^2)) / (L+sqrt(y+L^2)) )
 else PPTSmith <- 0
return(PPTSmith) 
}

PPT.Steele <- function(Ls,k,z,Pmax,Lopt){
 L <- Ls*exp(-k*z)
 if(Ls >0) PPTSteele <- ((Pmax*exp(1))/(k*z)) * ( exp(-L/Lopt) - exp(-Ls/Lopt) )
 else PPTSteele <- 0
 return(PPTSteele)
}

PPrates.depth <- function(z, param,sw.plot=T) {
# param light: k att coeff, Ls light subsurface
k <- param[1]; Ls <- param[2]
# param PP: Pmax,alpha,Lopt
Pmax <- param[3]; alpha <- param[4]; Lopt <- param[5]

# Light, Beer's law
L <- Ls*exp(-k*z)
# Smith & Steele model vs L
PPSmith <- PP.Smith(L,Pmax,alpha)
PPSteele <- PP.Steele(L,Pmax,Lopt)
PP <- cbind(PPSmith,PPSteele)
# depth average
PPTSmith <- PPT.Smith(Ls,k,z,Pmax,alpha)
PPTSteele <- PPT.Steele(Ls,k,z,Pmax,Lopt)
PPT <- cbind(PPTSmith,PPTSteele)
out <- data.frame(z,L,PP,PPT)

if(sw.plot==T){
 mat<- matrix(1:4,2,2,byrow=T)
 layout(mat,c(3.5,3.5),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
 # Light, Beer's law
 plot(L, -z, type="l",ylab="Depth [m]",xlab="Light [W/m2]")
 abline(h=0, col="grey")
 # Smith and Steele vs. L
 matplot(L,PP, type="l", col=1, ylab="P")
 legend("bottomright",legend=c("Smith","Steele"),lty=1:2)
 abline(h=0, col="grey")
 # smith and steele vs depth
 matplot(PP,-z, type="l", ylab="Depth [m]", xlab="PP",col=1)
 legend("topleft",legend=c("Smith","Steele"),lty=1:2)
 abline(h=0, col="grey")
 # depth avg
 matplot(PPT,-z,type="l", col=1,xlab="PP Depth Averaged",ylab="Depth [m]")
 legend("topleft",legend=c("Smith","Steele"),lty=1:2)
 abline(h=0, col="grey")
}
return(out) 
}


