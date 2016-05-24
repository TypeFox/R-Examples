CalcXP <-
function(p,v,z,zmin,zmax,novertex)
{


if( z[1]==z[2]){
   x=0
   return(x)
} 


if((v[1]==0)&&(v[2]==0)){
   x=Inf
   return(x)
} 

if(missing(zmin)){
	zmin = z[1]
	zmax = z[2]
}else{
	zmax=min(c(z[2],zmax))
	zmin=max(c(z[1],zmin))
}
if(missing(novertex)){novertex=0}

aa=SlopeInt(v,z)
g=aa[[1]]
v0=aa[[2]]

if( g==0){
   
   x=-p*v0*(zmin-zmax)/sqrt(1-p*p*v0*v0)
}else{
   

   if((z[1]<zmax)&&(zmax<z[2])&&(novertex==0)){
      
      
      
      
      
      
      v0=v[1]
      
      
      zs=zmax-z[1]
      v0gzs=v0+g*zs
      x=p*sqrt(g*zs*(v0+v0gzs))*v0gzs/g

   }else{

      
      
      
      p2g2=p*p*g*g
      v0rcpg=v0/g
      
      zplusv0g=zmin+v0rcpg
      xmin=sqrt(-p2g2*zplusv0g*zplusv0g+1+0i)
      
      zplusv0g=zmax+v0rcpg
      xmax=sqrt(-p2g2*zplusv0g*zplusv0g+1+0i)
      
      x=-(xmax-xmin)/(p*g)
          
      if(is.na(x)){
         stop('CalcXP: distance computation resulting in NaN!')
      } 
      
      
      
      
      
      
      
      if(Im(x)!=0){
         
         v0=v[1]
         zs=zmax-z[1]
         v0gzs=v0+g*zs
         x=p*sqrt(g*zs*(v0+v0gzs))*v0gzs/g
      }else{
      x=Re(x)	 
      }
   } 
} 

if(Im(x)!=0){
   x=NaN
} 

return(x)




}

