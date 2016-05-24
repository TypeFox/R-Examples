Kernel=function(u,type.Ker="Ker.norm"){
   tab=list("Ker.norm","Ker.cos","Ker.epa","Ker.tri","Ker.quar","Ker.unif")
   type.i=pmatch(type.Ker,tab)
   if (is.na(type.i))   Ker=dnorm(u)
     else {
     if (type.i==1)   Ker=dnorm(u)
     if (type.i==2)   Ker=ifelse(abs(u)<=1,pi/4*(cos(pi*u/2)),0)
     if (type.i==3)   Ker=ifelse(abs(u)<=1,0.75*(1-u^2),0)
     if (type.i==4)   Ker=ifelse(abs(u)<=1,35/32*(1-u^2)^3,0)
     if (type.i==5)   Ker=ifelse(abs(u)<=1,15/16*(1-u^2)^2,0)
     if (type.i==6)   Ker=ifelse(abs(u)<=1,0.5,0)
}
  return(Ker)
}

