CalcXPsum <-
function(p,v,z,zmin,zmax,novertex)
{
dist=NaN
segx=NaN
segz=NaN
res=NULL
resz=NULL



pepsilon=1e-10; 

anz=length(v)-1; 

done=0;
cnt=1;

while(!done){
   if (z[cnt]>=zmax){
      
      done=1
      res[cnt]=-1
      resz[cnt]=-1
   }else{
      if( (zmin>z[cnt])&(zmax>z[cnt])){
         res[cnt]=-1
         resz[cnt]=-1
      }else{
         
         layerv=c(v[cnt],v[cnt+1]) 
         layerz=c(z[cnt],z[cnt+1]) 
         if(novertex==1){
            res[cnt]=CalcXP(p,layerv,layerz,zmin,zmax,1)
            resz[cnt]=min(c(max(layerz),zmax))
         }else{
            res[cnt]=CalcXP(p,layerv,layerz,zmin,zmax)
            resz[cnt]=min(c(max(layerz),zmax))
         } 
      } 
   } 
   cnt=cnt+1 
} 

indies=which(res!=-1)
if(length(indies)){
   res=res[indies]
   resz=resz[indies]
}else{
   
   
   if(!is.null(res)){
      if(is.null(which(res!=-1))){
         res=c(0,0) 
         resz=c(0,0)
      } 
   } 
} 
indies=which(is.na(res))
if(length(indies)>0 ){
   indies=1:(indies[1]-1)
}else{
   indies=1:length(res)
} 



if( !length(indies)){
   
   dist=NaN
   segx=NaN
   segz=NaN
}else{
   dist=sum(res[indies]);  
   if( sum(res)==0){
      
      
      if( p<pepsilon){
         segx=c(0, res[indies]) 
         segz=c(zmin, resz[indies]); 
      }else{
          segx=res
          segz=res
      }
   }else{
      segx=c(0, res[indies]); 
      segz=c(zmin, resz[indies]) 
   } 
} 


if( max(segz)<zmax){
   
   dist=NaN
   segx=NaN
   segz=NaN
}
return(list(dist=dist,segx=segx,segz=segz))
}

