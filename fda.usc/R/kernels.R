Ker.norm=function(u){dnorm(u)}
Ker.cos=function(u){ifelse(abs(u)<=1,pi/4*(cos(pi*u/2)),0)}
Ker.epa=function(u){ifelse(abs(u)<=1,0.75*(1-u^2),0)}
Ker.tri=function(u){ifelse(abs(u)<=1,35/32*(1-u^2)^3,0)}
Ker.quar=function(u){ifelse(abs(u)<=1,15/16*(1-u^2)^2,0)}
Ker.unif=function(u){ifelse(abs(u)<=1,0.5,0)}

#skernel.norm<-function(x){sum(dnorm(x))}
rkernel<-function(x,y,Ker=Ker.norm){
             s1=sum(Ker(x)*y,na.rm=TRUE)
             ss=sum(Ker(x),na.rm=TRUE)
return(s1/ss)}

#distL2=function(x1,x2){mean((x1-x2)^2)} # Not used (see metri.lp)

AKer.norm=function(u)  {ifelse(u>=0,2*dnorm(u),0)}
AKer.cos=function(u)  {ifelse(u>=0,pi/2*(cos(pi*u/2)),0)}
AKer.epa=function(u) {ifelse(u>=0 & u<=1,1.5*(1-u^2),0)}
AKer.tri=function(u) {ifelse(u>=0 & u<=1,35/16*(1-u^2)^3,0)}
AKer.quar=function(u) {ifelse(u>=0 & u<=1,15/8*(1-u^2)^2,0)}
AKer.unif=function(u) {ifelse(u>=0 & u<=1,1,0)}



IKer.norm=function(u){
	kaux=function(u){integrate(Ker.norm,-1,u)$value}
	mapply(kaux,u)
}
#IKer.norm2=function(u){
#	kaux=pnorm(u)
#		mapply(kaux,u) falta una factor
#}

IKer.cos=function(u){
	kaux=function(u){integrate(Ker.cos,-1,u)$value}
	mapply(kaux,u)
}
#IKer.cos2=function(u){ifelse(abs(u)<=1,0.5*(sin(pi*u/2)),0)}


IKer.tri=function(u){
	kaux=function(u){integrate(Ker.tri,-1,u)$value}
	mapply(kaux,u)
}
IKer.quar=function(u){
	kaux=function(u){integrate(Ker.quar,-1,u)$value}
	mapply(kaux,u)
}
IKer.unif=function(u){
	kaux=function(u){integrate(Ker.unif,-1,u)$value}
	mapply(kaux,u)
}

IKer.epa=function(u){
	kaux=function(u){integrate(Ker.epa,-1,u)$value}
	mapply(kaux,u)
}

Kernel.integrate=function(u,Ker=Ker.norm,a=-1){
    kaux=function(u){integrate(Ker,a,u)$value}
    mapply(kaux,u)
}



