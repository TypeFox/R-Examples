# height vs diameter alometric relation
DH.allom <- function(D,param){
 H0 <- 1.37
 Hmax <- param[1]; b3 <- param[2]; b4 <- param[3]
 H <- H0 + (Hmax-H0)*(1 - exp(b3*D))^b4
 return(H)
}

Dl.allom <- function(D,param){
 c <- param[1]; fhc <- param[2]; d <- param[3]
 l <- c*(D)^d*fhc
 return(l)
}
  
dD.rate <- function(D,G,Dmax,param.H,param.l){
  H0=1.37
  H <- DH.allom(D,param.H)
  l <- Dl.allom(D,param.l)
  Hmax <- param.H[1]; b3 <- param.H[2]; b4 <- param.H[3]
  nume <- l*(1 - (D*H)/(Dmax*Hmax))
  phi <- 100*(2*H-D*(Hmax-H0)*b3*b4*exp(b3*D)*(1-exp(b3*D))^(b4-1))
  deno <- D*phi
  dD <- G*nume/deno
 return(dD)
} 

dD.dt.analysis <- function(D,G,Dmax,param.H,param.l,D0,t,pdfout=F){
H0=1.37
mat<- matrix(1:2,2,1,byrow=T)
layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

dD.dt <- matrix(nrow=length(D),ncol=length(G)) 
for(i in 1:length(G))
dD.dt[,i] <- dD.rate(D,G[i],Dmax,param.H,param.l)
matplot(D,dD.dt,type="l",col=1,xlab="Diameter D [cm]",ylab="Diameter increment dD/dt [cm/y]")
legend("topleft",legend=paste("Gmax=",G),lty=1:length(G),col=1)
mtext(side=1,line=-1,"Optimum conditions g(E(t))=1 for all t")

Dt <- matrix(nrow=length(t),ncol=length(G))
for (j in 1:length(G)){
Dt[1,j] <- D0  
for(i in 2:length(t))
Dt[i,j] <- Dt[i-1,j] + dD.rate(Dt[i-1,j],G[j],Dmax,param.H,param.l)
}
matplot(t,Dt,type="l",col=1,xlab="Time [yr]",ylab="Diameter D [cm]")
legend("bottomright",legend=paste("Gmax=",G),lty=1:length(G),col=1)
mtext(side=1,line=-1,"Optimum conditions g(E(t))=1 for all t")

if(pdfout==F) win.graph()
 mat<- matrix(1:2,2,1,byrow=T)
 layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

 Ht <- DH.allom(Dt,param.H)
 matplot(t,Ht,type="l",col=1,xlab="Time [yr]",ylab="Height H [m]")
 legend("bottomright",legend=paste("Gmax=",G),lty=1:length(G),col=1)
 mtext(side=1,line=-1,"Optimum conditions g(E(t))=1 for all t")

 lt <- Dl.allom(Dt,param.l)
 matplot(t,lt,type="l",col=1,xlab="Time [yr]",ylab="Leaf area l [m2]")
 legend("bottomright",legend=paste("Gmax=",G),lty=1:length(G),col=1)
 mtext(side=1,line=-1,"Optimum conditions g(E(t))=1 for all t")
 return(list(D.dD=data.frame(D,dD.dt),t.Xt =data.frame(t,Dt,Ht,lt)))
}

  parab1 <- function(E,Eopt,p){
  ne <- length(E)
  F <- 1.0-p*(Eopt-E)^2
   for(i in 1:ne) if(F[i]<0) F[i] <- 0
   for(i in 1:ne) if(F[i]>Eopt) F[i] <- 1
  return(F)
 }

 env.factors.analysis <- function(nsp,splab,FE,Eopt,param){
 # parabolic
 Elab <- c("L","T","W","N")
 mat<- matrix(1:4,2,2,byrow=T)
 layout(mat,c(3.5,3.5),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
 for(k in 1:4){
  E <- seq(0,Eopt[k],0.01*Eopt[k]); ne <- length(E)
  F <- matrix(nrow=ne,ncol=nsp)
  for(j in 1:nsp) F[,j] <- FE[[k]](E,Eopt[k],param[[j]][k])
  matplot(E,F,type="l",col=1,xlab=Elab[k],ylab=paste("Factor F(",Elab[k],")"))
  legend("bottomright",legend=splab, lty=1:nsp,col=1)
}
}

 ge.growth <- function(nsp,splab,FE,Eopt,param.E,Ef,param.Ef,Dmax,Gmax,param.H,param.l,D0,nt){
 E <- matrix(nrow=nt,ncol=4)
 G <- matrix(nrow=nt,ncol=nsp); D <- G
 F <- structure(1:(nt*nsp*4), dim=c(nt,nsp,4))
 D[1,1:nsp] <- D0
 E[1,] <- 1; F[1,,] <- 1; G[1,] <- 1 
 for(i in 2:nt){
   for(k in 1:4) {
     E[i,k] <- Ef[[k]](1,param.Ef[[k]])
     for(j in 1:nsp) F[i,j,k] <- FE[[k]](E[i,k],Eopt[k],param.E[[j]][k])
    } #k
    for(j in 1:nsp) {
     G[i,j] <- F[i,j,1]*F[i,j,2]*min(F[i,j,3],F[i,j,4])
     D[i,j] <- D[i-1,j] + dD.rate(D[i-1,j],Gmax[j]*G[i,j],Dmax[j],param.H[[j]],param.l[[j]])
    } #j
} #t

 mat<- matrix(1:2,2,1,byrow=T)
 layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
 matplot(t,G,type="s",xlab="Time [yr]",ylab="gE",col=1)
 legend("topleft", lty=1:nsp, legend=splab)
 matplot(t,D,type="l",xlab="Time [yr]",ylab="D [cm]",col=1)
 legend("topleft", lty=1:nsp, legend=splab)
 return(list(round(E,2),round(F,2),round(data.frame(t,E,G,D),2)))
}

dD.dt.calibra <- function(D,dD,Dmax,param.H,param.l){
mat<- matrix(1:2,2,1,byrow=T)
layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

fmax <- max(dD,na.rm=T)
Dopt <- D[which(dD==fmax)]
fopt <- dD.rate(Dopt,1,Dmax,param.H,param.l)

Gmax <- round(fmax/fopt,0)
matplot(D,dD, type="p",col=1,pch=1,xlab="Diameter D [cm]",ylab="Diameter increment dD/dt [cm/y]",ylim=c(0,2))
dD.dt <- array() 
dD.dt <- dD.rate(D,Gmax,Dmax,param.H,param.l)
lines(D,dD.dt,col=1)
legend("topleft",legend=paste("Gmax=",Gmax),lty=1,col=1)
mtext(side=1,line=-1,"Adjusted for Optimum conditions g(E)=1")
return(Gmax)
}

  expon.alf <- function(E,p){
  ne <- length(E)
  F <- p[1]*(1.0-exp(-p[2]*(E-p[3])))
   for(i in 1:ne) if(F[i]<0) F[i] <- 0
   for(i in 1:ne) if(F[i]>1) F[i] <- 1
  return(F)
 }

 gauss.smf <- function(x,p) {
 # calculates response to soil water conditions
 # WetDay: x2,  # DryDay: x3
 y <- x
 for(i in 1:2){
  if(x[i]<0) x[i] <- 1
  y[i]<- 1-exp(-(x[i]/p[i])^2.00/2.00)
 }
 # combine as limiting factor  
 smf<- min(y[1],y[2])
 return(c(y,smf))
}





