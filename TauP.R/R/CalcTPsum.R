CalcTPsum <-
function(p,v,z,zmin,zmax,novertex){
ttime=NaN
res= 0 


anz=length(v)-1


done=0
cnt=1
while(!done){
   if( z[cnt]>=zmax){
      
      done=1
   }else{
      if((zmin>z[cnt])&(zmax>z[cnt])){
         res[cnt]=-1
      }else{
         
         layerv=c(v[cnt],v[cnt+1]) 
         layerz=c(z[cnt], z[cnt+1]) 

            res[cnt]=CalcTP(p,layerv,layerz,zmin,zmax,novertex)

      } 
   } 
   cnt=cnt+1 
} 

if( sum(is.na(res))>0){
   res=NaN
   return(ttime)
} 

indies=which(Re(res)>-1) 
if(length(indies)){
   res=res[indies]
}else{
   
   
   if(length(res)){
      if(!length(which(res!=-1))){
         res=0 
      } 
   } 
} 
indies=which(is.na(res))
if(length(indies)){
   indies=1:(indies[1]-1)
}else{
   indies=1:length(res)
} 

ttime=sum(res[indies])  

return(ttime)
}

