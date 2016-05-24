CalcTP <-
function(p,v,z,zmin,zmax,novertex=0){

### Define internal function:
EvalT=function(z,g,p2g,v0rcpg)
{
zplusv0g=z+v0rcpg
denominator=sqrt(-p2g*g*zplusv0g*zplusv0g+1+0i)
if (Re(denominator)==0){
   t=0 
   return(t)
}
t=atanh(1/denominator+0i) 
return(t)
}
### Done defining internal function
  
if( z[1]==z[2]){
   t=0
   return(t)
} 
if( (v[1]==0)&(v[2]==0)){
   t=Inf
   return(t)
}
if(missing(zmin)){
   zmin=z[1]
   zmax=z[2]
}else{
   zmax=min(c(z[2],zmax))
   zmin=max(c(z[1],zmin))
 } 
aa=SlopeInt(v,z)
g=aa[[1]]
v0=aa[[2]]
if( g==0){
   
   t=(zmax-zmin)/sqrt(v0*(1-p*p*v0)) 
}else{
   
   if((z[1]<zmax)&(zmax<z[2])&!novertex){
      
      
      
      
      
      v0=v[1]
      
      pv0=p*v0
      
      t=log((1+sqrt(1-pv0*pv0))/pv0)/g
   }else{
      
      
      p2g=p*p*g
      v0rcpg=v0/g
      
      tmin=EvalT(zmin,g,p2g,v0rcpg)
      
      tmax=EvalT(zmax,g,p2g,v0rcpg)
      
      t=(-1/g)*(tmax-tmin)
      
      if( is.infinite(Re(t))){
         t=Inf 
      } 
      
      
      if( Im(t)!=0){
         
         
         v0=v[1]
         pv0=p*v0
         t=log((1+sqrt(1-pv0*pv0))/(pv0))/g
      } 
   } 
} 
if(Im(t)!=0){
   t=NaN
} 
return(Re(t))
}

