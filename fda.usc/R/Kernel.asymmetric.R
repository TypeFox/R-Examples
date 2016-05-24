Kernel.asymmetric=function(u,type.Ker="AKer.norm"){
   tab=list("AKer.norm","AKer.cos","AKer.epa","AKer.tri","AKer.quar","AKer.unif")
   type.i=pmatch(type.Ker,tab)
   if (is.na(type.i))   Ker=ifelse(u>=0,2*dnorm(u),0)
     else {
     if (type.i==1)   Ker=ifelse(u>=0,2*dnorm(u),0)
		 if (type.i==2)   Ker=ifelse(u>=0,pi/2*(cos(pi*u/2)),0)
     if (type.i==3)   Ker=ifelse(u>=0 & u<=1,1.5*(1-u^2),0)
     if (type.i==4)   Ker=ifelse(u>=0 & u<=1,35/16*(1-u^2)^3,0)
     if (type.i==5)   Ker=ifelse(u>=0 & u<=1,15/8*(1-u^2)^2,0)
     if (type.i==6)   Ker=ifelse(u>=0 & u<=1,1,0)
}
  return(Ker)
}

