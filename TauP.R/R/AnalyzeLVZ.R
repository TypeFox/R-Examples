AnalyzeLVZ <-
function(v,vsec,z,rp)
{

zepsilon=0.001; 

criticalz=NULL 
newz=NULL 
newv=NULL 
newvsec=NULL 


zerolist=which(v==0)

if(!is.null(zerolist)){
    
    
    
    v[zerolist]=vsec[zerolist]}




a=TransformS2Fz(v,z,rp) 
vf=a[[1]]
zf=a[[2]]
vsecf=TransformS2Fz(vsec,z,rp)[[1]] 





sampanz=length(z)
inlvz=0 
lastabove=NULL 
firstbelow=NULL 
vtop=NULL 
vsectop=NULL 


for(sampcnt in 2:sampanz){ 
    if((inlvz==0)&(vf[sampcnt]<vf[sampcnt-1])){
       
       inlvz=1 
       lastabove=c(lastabove, sampcnt-1) 
       vtop=c(vtop,vf[sampcnt-1]) 
       vsectop=c(vsectop,vsecf[sampcnt-1])
       }
    
    if((inlvz==1)& prod(c(TRUE,(vf[sampcnt]>=vtop[length(vtop)])))){
       
       inlvz=0 
       firstbelow=c(firstbelow, sampcnt)
       }
    }



bottomanz=length(firstbelow) 
bottomz=firstbelow+NaN 
bottomv=bottomz 
bottomvsec=bottomz 
for(lvzcnt in 1:bottomanz){
    
    vbelow=vf[firstbelow[lvzcnt]] 
    vabove=vf[firstbelow[lvzcnt]-1]; 
    vsecbelow=vsecf[firstbelow[lvzcnt]] 
    vsecabove=vsecf[firstbelow[lvzcnt]-1]; 
    zbelow=zf[firstbelow[lvzcnt]] 
    zabove=c(0,zf)[firstbelow[lvzcnt]]; 

    if(zabove!=zbelow){
        
        
        
        bottomz[lvzcnt]=approx(c(vabove, vbelow),c(zabove, zbelow),vtop[lvzcnt],'linear')$y
        bottomv[lvzcnt]=vtop[lvzcnt]
	}
       	
        
        
   
    

    
    
    vsecbelow=vsecf[firstbelow[lvzcnt]] 
    vsecabove=vsecf[firstbelow[lvzcnt]-1] 

    if(zabove!=zbelow){
        
        
 
        bottomvsec[lvzcnt]=approx(c(zabove,zbelow),c(vsecabove,vsecbelow),bottomz[lvzcnt],'linear')$y
	}
}

remain=which(!is.na(bottomz))
bottomz=bottomz[remain]
bottomv=bottomv[remain]
bottomvsec=bottomvsec[remain]

topz=zf[lastabove]

newz=bottomz 
      
      
newv=bottomv 
      
      
newvsec=bottomvsec
             
 
topz=TransformF2Sz(topz,topz,rp)[[2]]
bottomz=TransformF2Sz(bottomz,bottomz,rp)[[2]]
dmy=newz; 
a=TransformF2Sz(newv,newz,rp)
newv=a[[1]]
newz=a[[2]]
newvsec=TransformF2Sz(newvsec,dmy,rp)[[1]]



bottomz=bottomz+zepsilon

criticalz=c(topz,bottomz)

return(list(newv=newv,newvsec=newvsec,newz=newz,criticalz=criticalz))
}

