ConvP2Vdepth <-
function(p,v,r,h,rp,discons)
{

  # Define interpolation function:


  Interp2pt=function(x,y,xnew)
{


sorter=order(x)
x=sort(x)
y=y[sorter]
infties=which(is.infinite(x))

    if(length(infties)==0){
        
        if( x[2]-x[1]!=0){
            
            ynew=((x[2]-xnew)*y[1]+(xnew-x[1])*y[2])/(x[2]-x[1])
        }else{
            
            ynew=NaN
        } 
    }else if(length(infties)==1){
        
        if( is.infinite(x[1])){
           
           ynew=y[2]
        }else{
           
           ynew=y[1]
        } 
    }else if(length(infties)==2){
        
        
        ynew=NaN
    } 

return(ynew)
}

### Done defining interpolation function

sloepsilon=1e-6

zprecision=1e-6 

p = p * 180/pi


res=NaN



aa=TransformS2Fz(v,rp-r,rp)
vflat=aa[[1]]
zflat=aa[[2]]

pflat=p/rp 
hflat=TransformS2Fz(rp-h,rp-h,rp)[[2]]


slo=1/vflat


sourceslo=LinInterp(zflat,slo,hflat,'all'); 


    if(length(sourceslo)==1){
         
         
         if( abs(pflat-sourceslo)<sloepsilon){
            rpd=h
            return(rpd)
         } 
    }else if(length(sourceslo)==2){
         
         
         
         if( (pflat<=sourceslo[1])&(pflat>=sourceslo[2])){ 
							   
							   
            rpd=h
            return(rpd)
         } 
    }else{
         
         
         stop('ConvP2Vdepth: unexpected number of slownesses at source depth')
    } 



sampleanz=length(slo) 
topz=zflat[1:(sampleanz-1)] 
botz=zflat[2:sampleanz]  
topslo=slo[1:(sampleanz-1)]  
botslo=slo[2:sampleanz]   

	
candidatelayers=which(((topslo<=pflat)&(pflat<=botslo))|((topslo>=pflat)&(pflat>=botslo)))
if(length(candidatelayers)){
   
                 
                 
    
    
    
    
    topcandidate=topz[candidatelayers] 
    botcandidate=botz[candidatelayers] 
    topslocandidate=topslo[candidatelayers] 
    botslocandidate=botslo[candidatelayers] 

    
    deepenough=which((botcandidate>hflat) & 
                    (topslocandidate>=botslocandidate)) 
                
    if(length(deepenough)){
        
        

        topcandidate=topcandidate[deepenough]
        botcandidate=botcandidate[deepenough]
        topslocandidate=topslocandidate[deepenough]
        botslocandidate=botslocandidate[deepenough]


        
        
        topcandidate=topcandidate[1]
        botcandidate=botcandidate[1]
        topslocandidate=topslocandidate[1]
        botslocandidate=botslocandidate[1]

        
        
        turndepthflat=Interp2pt(c(topslocandidate, botslocandidate),c(topcandidate, botcandidate),pflat)

        
        res=TransformF2Sz(turndepthflat,turndepthflat,rp)[[2]]
        res=rp-res 

    }else{
        
        res=NaN
    }

}else{
   
   res=NaN
} 

if(!is.na(res)){
   
   dr=abs(discons-res)
   
   indy=which(dr<=zprecision)
   if(length(indy)){
      
      
      res=discons[indy]
   } 
} 



return(res)
}

