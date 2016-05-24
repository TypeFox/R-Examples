FindPrange <-function(phase,imodel,h,dangle){
  angles=NULL
  radian=pi/180
  
  phasev=strsplit(phase,'')[[1]] 

  if(!is.na(imodel$cmb) && h > imodel$cmb){stop('Focal depth must be above the CMB')}
  if(!is.na(imodel$cmb) && h == imodel$cmb && phasev[2] %in% c('c', 'K')){stop("Focal depth must be strictly above CMB for this phase")}
  
  
if( phasev[1] %in% c('p','s') ) {
   
   
  if(h == 0){
    stop('Depth phases do not exist for h = 0')
  }
   depthleg=phasev[1] 
   surfaceleg=phasev[2:length(phasev)] 
   
   
   
   
   aa=FindPrange(paste(surfaceleg,collapse=''),imodel,0,10)
   angles=aa[[1]]
   minangle=aa[[2]]
   maxangle=aa[[3]]
   
   
   
   raypminangle=ConvAng2p(paste(surfaceleg, collapse = ''),0,minangle,model = NULL, imodel$vp[1],imodel$vs[1],imodel$rp)
   raypmaxangle=ConvAng2p(paste(surfaceleg, collapse = ''),0,maxangle,model = NULL, imodel$vp[1],imodel$vs[1],imodel$rp)
   
   
   focus=InterpModel(imodel,h,'simple')


   minangle=ConvP2Ang(depthleg,h,raypminangle,vp=focus$vp[1],vs=focus$vs[1],rp=focus$rp)
   maxangle=ConvP2Ang(depthleg,h,raypmaxangle,vp=focus$vp[1],vs=focus$vs[1],rp=focus$rp)
   



   
   
   
   
   
   
   
   minangle=minangle[!(minangle<90)]
   maxangle=maxangle[!(maxangle<90)]
   
   
   
   
   if( sum(is.na(minangle)) ){minangle=90}
   if( sum(is.na(maxangle)) ){maxangle=90}
   
   
   
   
   if( minangle>maxangle ){
      swapper=minangle
      minangle=maxangle
      maxangle=swapper
   } 
   
   
   

   if( is.na(minangle+dangle+maxangle) ){
      angles=NULL
   }else{
      angles=seq(from = minangle, by = dangle, to = maxangle)
      if( angles[length(angles)]!=maxangle ){angles=c(angles,maxangle)}
   } 

   
   
   
   return(list(angles=angles,minangle=minangle,maxangle=maxangle))
   
   
} 





if(!('criticalrays' %in% names(imodel)) ){
    
    imodel=ImproveModel(imodel)[[1]]
} 

imodel$criticalrays$z = round(imodel$criticalrays$z, 3)


focus=InterpModel(imodel,h,'simple')


if( prod(focus$vs == 0) ){ #5/16/11: changed from focus$vs == 0 to prod(focus$vs) == 0 because it produced warnings on discontinuities.  used prod instead of sum to permit sources on the CMB but not inside the outer core. 
   angles=NULL
   return(list(angles=angles,minangle=minangle,maxangle=maxangle))
} 


psurf=imodel$rp/imodel$vp[1]
ssurf=imodel$rp/imodel$vs[1]
psurf=psurf * pi/180
ssurf=ssurf * pi/180
psurfangle=ConvP2Ang('P',0,psurf,vp=imodel$vp[1],vs=imodel$vs[1],rp=imodel$rp)
ssurfangle=ConvP2Ang('S',0,ssurf,vp=imodel$vp[1],vs=imodel$vs[1],rp=imodel$rp)




psteepangle=0.1 
ssteepangle=0.1 


aa=ConvVdepth2p(imodel,imodel$z[length(imodel$z)-1])
psteep=aa[[1]]
ssteep=aa[[2]]
psteepangle=min(ConvP2Ang('P',h,psteep,vp=focus$vp[1],vs=focus$vs[1],rp=imodel$rp));
ssteepangle=min(ConvP2Ang('S',h,ssteep,vp=focus$vp[1],vs=focus$vs[1],rp=imodel$rp));




if( !is.na(imodel$cmb) ){
   
   
   
   
   indies=which(0 == round(imodel$criticalrays$z - imodel$cmb, 3))
   pcmb=imodel$criticalrays$p[indies] 
   scmb=imodel$criticalrays$s[indies] 
   
   
   
   
   indies=which(!is.na(pcmb))
   pcmb=pcmb[indies[1]]
   indies=which(!is.na(scmb)) 
   scmb=scmb[indies[1]]



   
   
   
   
   
   if( h<=imodel$cmb ){
       
       vp=focus$vp[length(focus$vp)]
       vs=focus$vs[length(focus$vs)]
   }else{
       
       vp=focus$vp[1]
       vs=focus$vs[1]
   } 
   
   
   pcmbangle=ConvP2Ang('P',h,pcmb,vp=vp,vs=vs,rp=imodel$rp);
   scmbangle=ConvP2Ang('S',h,scmb,vp=vp,vs=vs,rp=imodel$rp);
   
   
   if( h<=imodel$cmb ){
       
       pcmbangle=pcmbangle[length(pcmbangle)] 
       scmbangle=scmbangle[length(scmbangle)] 
   }else{
       
       pcmbangle=pcmbangle[1] 
       scmbangle=scmbangle[1] 
   } 
   
   
}else{
   
   
   pcmb=psteep
   scmb=ssteep
   pcmbangle=psteepangle
   scmbangle=ssteepangle
} 


if( !is.na(imodel$icb) ){
   
   
   
   
   
   indies=which(0 == round(imodel$criticalrays$z - imodel$icb, 3))
   picb=imodel$criticalrays$p[indies] 
   sicb=imodel$criticalrays$s[indies] 
   
   
   
   
   indies=which(!is.na(picb))
   picb=picb[indies[1]]
   indies=which(!is.na(sicb)) 
   sicb=sicb[indies[1]]
   

   
   
   
   
   
   if( h<=imodel$icb ){
       
       vp=focus$vp[length(focus$vp)]
       vs=focus$vs[length(focus$vs)]
   }else{
       
       vp=focus$vp[1]
       vs=focus$vs[1]
   } 
   
   
   picbangle=ConvP2Ang('P',h,picb,vp=vp,vs=vs,rp=imodel$rp)
   sicbangle=ConvP2Ang('S',h,sicb,vp=vp,vs=vs,rp=imodel$rp)
   skicbangle=ConvP2Ang('S',h,picb,vp=vp,vs=vs,rp=imodel$rp) 
   
   
   if( h<=imodel$cmb ){
       
       picbangle=picbangle[length(picbangle)] 
       sicbangle=sicbangle[length(sicbangle)] 
       skicbangle=skicbangle[length(skicbangle)]
   }else{
       
       picbangle=picbangle[1] 
       sicbangle=sicbangle[1] 
       skicbangle=skicbangle[1]
   } 

}else{
   
   
   picb=psteep
   sicb=ssteep
   picbangle=psteepangle
   sicbangle=ssteepangle
   skicbangle=psteepangle
} 





if( !is.na(imodel$cmb) ){
    
    
    cmbindies=which(imodel$z==imodel$cmb)
  
    
    
    
    vpcmb=imodel$vp[cmbindies]
    
    
    
    
    
    
    cmbreflectp=(imodel$rp-imodel$cmb)/vpcmb[2]
    
    
    cmbreflectp=cmbreflectp * pi/180
    
    
    
    pcmbreflectangle=ConvP2Ang('P',h,cmbreflectp,vp=vp,vs=vs,rp=imodel$rp)
    scmbreflectangle=ConvP2Ang('S',h,cmbreflectp,vp=vp,vs=vs,rp=imodel$rp)
    
    
    if( h<=imodel$cmb ){
        
        pcmbreflectangle=pcmbreflectangle[length(pcmbreflectangle)] 
        scmbreflectangle=scmbreflectangle[length(scmbreflectangle)]; 
    }else{
        
        pcmbreflectangle=pcmbreflectangle[1] 
        scmbreflectangle=scmbreflectangle[1] 
    } 
    
}else{
    
    scmbreflectangle=NaN
    pcmbreflectangle=NaN
    cmbreflectp=NaN
} 



if( !is.na(imodel$icb) ){
    
    
    icbindies=which(imodel$z==imodel$icb)
    
    
    
    
    vpicb=imodel$vp[icbindies[2]]
    vsicb=imodel$vs[icbindies[2]]
    
    
    picbshallowp=(imodel$rp-imodel$icb)/vpicb
    sicbshallowp=(imodel$rp-imodel$icb)/vsicb
    
    
    picbshallowp=picbshallowp * pi/180
    sicbshallowp=sicbshallowp * pi/180
    
    
    picbshallowangle=ConvP2Ang('P',h,picbshallowp,vp=vp,vs=vs,rp=imodel$rp)
    sicbshallowangle=ConvP2Ang('S',h,sicbshallowp,vp=vp,vs=vs,rp=imodel$rp)
    skicbshallowangle=ConvP2Ang('S',h,picbshallowp,vp=vp,vs=vs,rp=imodel$rp) 
    
    
    if( h<=imodel$icb ){
        
        picbshallowangle=picbshallowangle[length(picbshallowangle)]; 
        sicbshallowangle=sicbshallowangle[length(sicbshallowangle)] 
        skicbshallowangle=skicbshallowangle[length(skicbshallowangle)]
    }else{
        
        picbshallowangle=picbshallowangle[1] 
        sicbshallowangle=sicbshallowangle[1] 
        skicbshallowangle=skicbshallowangle[1]
    } 

}else{
    
    sicbshallowangle=NaN
    picbshallowangle=NaN
    skicbshallowangle=NaN
    skicbshallowp=NaN
    picbshallowp=NaN
    sicbshallowp=NaN
} 







aa=StripRepetitions(phase)
phase=aa[[1]]
repetitions=aa[[2]]
    

   if(phase == 'P'){        
        minangle=pcmbangle
        if(is.na(minangle)){minangle=psteepangle}
        if(h>0){
	  maxangle=180
        }else{
          maxangle=90
        } 
   }else if(phase %in% c('PP','PPP')){
        minangle=pcmbangle
        if(is.na(minangle)){
           minangle=psteepangle
        } 
        maxangle=90  
   }else if(phase %in% c('PcP','PcPPcP')){
        minangle=psteepangle
        maxangle=pcmbangle
   }else if(phase %in% c('PKP','PKKP','PKPPKP')){
        minangle=psteepangle
        maxangle=min(c(pcmbreflectangle,pcmbangle))
   }else if(phase == 'PKiKP'){
        minangle=psteepangle
        maxangle=picbangle
   }else if(phase %in% c('PKIKP','PKIKPPKIKP')){
        minangle=psteepangle
        maxangle=picbshallowangle
   }else if(phase == 'S'){
        minangle=scmbangle
        if(is.na(minangle)){
           minangle=ssteepangle
        } 
        if( h>0){
          maxangle=180
        }else{
          maxangle=90
        } 
   }else if(phase %in% c('SS','SSS')){
        minangle=scmbangle
        if( is.na(minangle)){
           
           minangle=ssteepangle
        } 
        maxangle=90 
   }else if(phase %in% c('ScS','ScSScS')){
        minangle=ssteepangle
        if(is.na(imodel$icb)){
           
           
           
           minangle=psteepangle
        }else{
           minangle=ssteepangle
        } 
        maxangle=scmbangle
   }else if(phase %in% c('SKS','SKKS','SKSSKS')){      
        if( imodel$vs[length(imodel$vs)-1]==0){
            
            
            
            minangle=psteepangle
        }else{
            
            minangle=ssteepangle
        } 
        maxangle=min(c(scmbreflectangle,scmbangle))
   }else if(phase == 'SKiKS'){        
        minangle=ssteepangle
        maxangle=skicbangle 
   }else if(phase == 'SKIKS'){
        minangle=ssteepangle
        maxangle=skicbshallowangle 
   }else if(phase == 'PS'){
        if(!is.na(imodel$cmb)){
            
            minangle=ConvP2Ang('P',h,scmb,vp=vp,vs=vs,rp=imodel$rp)
        }else{
            
            
            minangle=ConvP2Ang('P',h,ssteep,vp=focus$vp[length(focus$vp)],vs=focus$vs[length(focus$vs)],rp=imodel$rp[length(imodel$rp)])
        } 
        minangle=minangle[length(minangle)]
        if( is.na(minangle)){
           minangle=psteepangle
        } 
        maxangle=min(c(psurfangle, 90)) 
   }else if(phase == 'SP'){
        minangle=scmbangle
        if(is.na(minangle)){
           minangle=ssteepangle
        } 
        maxangle=min(c(psurfangle, 90)) 
   }else if(phase %in% c('PcS','PcSPcS')){
        minangle=psteepangle 
        maxangle=pcmbangle
   }else if(phase %in% c('ScP','ScPScP')){
        minangle=psteepangle 
        maxangle=scmbangle
   }else{
       print(paste('MKSMARTTAKEOFF: unknown phase ', phase, ', using default.'))
       minangle=psteepangle
       if( h>0){
          maxangle=180
       }else{
          maxangle=90
       } 
   } 




if( is.na(minangle+dangle+maxangle)){
   angles=NULL
}else{
   angles=seq(from=minangle,by=dangle,to=maxangle)
   if( angles[length(angles)]!=maxangle){
      angles=c(angles, maxangle)
   } 
} 


return(list(angles=angles,minangle=minangle,maxangle=maxangle))
}

