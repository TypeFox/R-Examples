FindTime4p <-
function(phase,h,p,model,anglemode='rayparm',takeoff=NULL)
{
tt=NaN
vdep=NaN
resp=NaN
radian=pi/180
zprecision=1e-6 
angleprecision=1e-6
resp=NULL
if(anglemode=='rayparm'){
  resp=p
}else if(anglemode=='angle'){
  if(is.null(takeoff)){takeoff = p}
   p=ConvAng2p(phase,h,takeoff,model)
   resp=p
}else{
   stop("FindTime4p: anglemode must be one of 'rayparm' or 'angle'")
}
if(is.null(takeoff)){takeoff=-Inf}
imagepsilon=1e-6 
if(abs(Im(takeoff))>imagepsilon){
   
   tt=NaN
   return(list(tt = tt, vdep = vdep, resp = resp))
}else{
   
   takeoff=Re(takeoff)
} 
if( (is.infinite(p))|(is.na(p))|(is.na(takeoff))){
    dist=NaN
    return(list(tt = tt, vdep = vdep, resp = resp))
} 
modelsav=model
psav=p
hsav=h
firstchar=substr(phase,1,1)
if(firstchar %in% c('p','s')){
   
   remphase=substr(phase,2,nchar(phase)) 
   
   
   
   if( h==0){
      
      stop('MKTIM4P: Depth phases do not exist for surface foci!')
   } 
   if( (takeoff<90)&(!is.infinite(takeoff))){
      return(list(tt = tt, vdep = vdep, resp = resp))
   } 
   
   
   if( is.infinite(takeoff)){
      dlp=ConvP2Ang(toupper(firstchar),h,p,model)
      dlp=dlp[1] 
      anglemode='angle'
   }else{
      dlp=takeoff
      anglemode='angle'
   } 
   aa=FindTime4p(toupper(firstchar),h,dlp,model,anglemode)
   dltt=aa[[1]]
   dlvdep=aa[[2]]
   dlresp=aa[[3]]
   
   
   if( is.infinite(takeoff)){
      
      
      remp=p
      anglemode='rayparm'
   }else{
      
      
      
      remp=ConvAng2p(firstchar,h,takeoff,model) 
      anglemode='rayparm'
   } 
   aa=FindTime4p(remphase,0,remp,model,anglemode)
   remtt=aa[[1]]
   remvdep=aa[[2]]
   remresp=aa[[3]]
   
   
   tt=dltt+remtt
   vdep=max(dlvdep,remvdep)
   resp=dlresp
   
   
   return(list(tt = tt, vdep = vdep, resp = resp))
}
aa=StripRepetitions(phase)
remphase=aa[[1]]
repetitions=aa[[2]]
if( repetitions>1){
   
   if( takeoff<=90){
       
       
       if( is.infinite(takeoff)){
          
          
          basep=p
          anglemode='rayparm'
       }else{
          
          
          
          basep=ConvAng2p(firstchar,h,takeoff,model) 
          anglemode='rayparm'
       } 
   aa=FindTime4p(remphase,h,basep,model,anglemode)
   basett=aa[[1]]
   basevdep=aa[[2]]
   baseresp=aa[[3]]
       
       if( is.infinite(takeoff)){
          
          
          remp=p
          anglemode='rayparm'
       }else{
          
          
          
          remp=ConvAng2p(firstchar,h,takeoff,model) 
          anglemode='rayparm'
       } 
       aa=FindTime4p(remphase,0,remp,model,anglemode)
   remtt=aa[[1]]
   remvdep=aa[[2]]
   remresp=aa[[3]]
       
       tt=basett + (repetitions - 1) * remtt
       resp=baseresp
       vdep=min(basevdep,remvdep) 
       
       return(list(tt = tt, vdep = vdep, resp = resp))
   
   
   }else{
       
       
       segx=NaN
       segy=NaN
       segtyp=NaN
       dist=NaN
       return(list(tt = tt, vdep = vdep, resp = resp))
       
   } 
       
} 
if( max(model$z)<model$rp){
   anz=length(model$z) 
   model$z=c(model$z, model$rp)
   model$vp=c(model$vp,model$vp[anz])
   model$vs=c(model$vs, model$vs[anz])
   model$rho=c(model$rho, model$rho[anz])
   model$qp=c(model$qp,model$qp[anz])
   model$qs=c(model$qs,model$qs[anz])
} 
vp=model$vp
vs=model$vs
rp=model$rp
r=model$rp-model$z
h=rp-h
aa=TransformS2Fz(vp,rp-r,rp)
vpflat=aa[[1]]
zflat=aa[[2]]
vsflat=TransformS2Fz(vs,rp-r,rp)[[1]]
cmbspher=model$cmb 
icbspher=model$icb 
alld=c(model$conr, model$moho, model$d410,model$d520, model$d660, model$cmb, model$icb, model$dz)
alld=TransformS2Fz(alld,alld,model$rp)[[2]]
model$conr=alld[1]
model$moho=alld[2]
model$d410=alld[3]
model$d520=alld[4]
model$d660=alld[5]
model$cmb=alld[6]
model$icb=alld[7]
model$dz=alld[8:length(alld)]
disconradii=FindDiscon(modelsav)
lowermostflat=zflat[length(zflat)-1]
if( is.na(model$cmb)){
   model$cmb=lowermostflat
   modelsav$cmb=modelsav$z[length(modelsav$z)-1]
   cmbspher=modelsav$z[length(modelsav$z)-1]
} 
if( is.na(model$icb)){
   model$icb=lowermostflat
   modelsav$icb=modelsav$z[length(modelsav$z)-1]
   icbspher=modelsav$z[length(modelsav$z)-1]
} 
  
if(phase=='PS'){
  aa=FindTime4p('P',hsav,psav,modelsav,takeoff=takeoff)
  ptt=aa[[1]]
  pvdep=aa[[2]]
  if( takeoff<=90){
    reflecttakeoff=ConvP2Ang('S',0,psav,modelsav)
    aa=FindTime4p('S',0,reflecttakeoff,modelsav,'angle') 
    stt=aa[[1]]
    svdep=aa[[2]]
    if( !is.na(stt)){
      if( (abs(svdep-modelsav$cmb)>zprecision)&(abs(pvdep-modelsav$cmb)>zprecision)){
        tt=ptt+stt
        vdep=min(svdep,pvdep)
      } 
    }else{
      tt=NaN
    } 
  }else{
    tt=NaN
  } 
}else if(phase == 'SP'){
  aa=FindTime4p('S',hsav,psav,modelsav,takeoff=takeoff)
  stt=aa[[1]]
  svdep=aa[[2]]
  if( takeoff<=90){
    reflecttakeoff=ConvP2Ang('P',0,psav,modelsav) 
    aa=FindTime4p('P',0,reflecttakeoff,modelsav,'angle') 
    ptt=aa[[1]]
    pvdep=aa[[2]]
    if( !is.na(ptt)){
      if( (abs(svdep-modelsav$cmb)>zprecision)&(abs(pvdep-modelsav$cmb)>zprecision)){
        tt=ptt+stt
        vdep=min(svdep,pvdep)
      } 
    }else{
      tt=NaN
    }
  }else{
    tt=NaN
  } 
  
  
  
}else if(phase %in% c('P','PP','PPP','PKP','PKPPKP','PKIKP','PKIKPPKIKP')){
  vdep=ConvP2Vdepth(p,vp,r,h,model$rp,disconradii) 
  if( is.na(vdep)&(takeoff<90) ) {
    tt=NaN
    return(list(tt = tt, vdep = vdep, resp = resp))
  }else{
    p=TransformS2Fp(p,rp)
    hflat=TransformS2Fz(c(1,1),rp-h,rp)[[2]]
    if( takeoff<=90){
      indies=which(vdep<=h) 
      vdep=max(vdep[indies]) 
      if( vdep<r[length(r)-1]){
        
        tt=Inf
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      if( is.null(vdep)){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      } 
      
      if(phase %in% c('PKP','PKPPKP')){
        if( is.na(model$cmb)){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
        if( vdep>=rp-cmbspher){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
        
      }else if(phase %in% c('PKIKP','PKIKPPKIKP')){
        if( is.na(model$cmb)){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
        if( is.na(model$icb)){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
        if( vdep>=rp-icbspher ){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        }
      }else if(phase %in% c('P','PP','PPP')){
        if( vdep-(rp-cmbspher)<zprecision){
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
      }
    } 
    tt2=CalcTPsum(p,vpflat,zflat,0,hflat,1) 
    if( takeoff<=90 ){
      vdepflat=TransformS2Fz(c(1,1),rp-vdep,rp)[[2]]
      tt1=CalcTPsum(p,vpflat,zflat,0,vdepflat,0)
      
      tt=2*tt1-tt2
      
      if(phase %in%  c('PP','PKPPKP','PKIKPPKIKP')){
        tt=tt+2*tt1
      }else if(phase == 'PPP'){
        tt=tt+4*tt1
      } 
    }else{
      if(phase=='P'){
        tt=tt2
      }else{
        tt=NaN
      } 
    } 
  } 
  
  
    
}else if(phase=='PKKP'){
  if( takeoff<=90){
    vdep=ConvP2Vdepth(p,vp,r,h,model$rp,disconradii) 
    if( is.na(vdep)){
      
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    }else{
      if( is.na(model$cmb)){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      indies=which(vdep<=h)
      vdep=max(vdep[indies]) 
      if( vdep<r[length(r)-1]){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      if( is.null(vdep)){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      } 
      p=TransformS2Fp(p,rp)
      aa=TransformS2Fz(c(1, 1),rp-c(vdep, h),rp)
      vdepflat=aa[[2]][1]
      hflat=aa[[2]][2]
      if( (vdepflat<=model$cmb)){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      } 
      
      tt1=CalcTPsum(p,vpflat,zflat,0,model$cmb,0)
      
      tt2=CalcTPsum(p,vpflat,zflat,0,hflat,1)
      
      ttk=CalcTPsum(p,vpflat,zflat,model$cmb,vdepflat,0)
      if(phase=='PKKP'){
        tt=2*tt1-tt2+4*ttk
      } 
    } 
  }else{
    tt=NaN
  } 
  
}else if(phase %in% c('PcP','PcPPcP')){
  if( takeoff<=90){
    
    maxp=ConvP2Vdepthinv(model$rp-cmbspher,model$vp,model$rp-model$z)
    if( p>maxp){
      
      tt=Inf
      return(list(tt = tt, vdep = vdep, resp = resp))
    }
    vdep=ConvP2Vdepth(p,vp,r,h,model$rp,disconradii) 
    indies=which(vdep<=h) 
    vdep=max(vdep[indies]) 
    if( is.na(model$cmb)){
      
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    }
    if( is.null(vdep)){
      
      vdep=TransformF2Sz(model$cmb+1,model$cmb+1,model$rp) 
    } 
    p=TransformS2Fp(p,rp)
    aa=TransformS2Fz(c(1,1),rp-c(vdep,h),rp)[[2]]
    vdepflat=aa[1]
    hflat=aa[2]
    if( vdepflat<model$cmb){
      
      tt=Inf
      return(list(tt = tt, vdep = vdep, resp = resp))
    } 
    tt1=CalcTPsum(p,vpflat,zflat,0,model$cmb,0)
    tt2=CalcTPsum(p,vpflat,zflat,0,hflat,1) 
    
    if(phase =='PcP'){
      tt=2*tt1-tt2}
    if(phase == 'PcPPcP'){
      tt=4*tt1-tt2}
  }else{
    tt=NaN
  } 
  
  
  
  
  
}else if(phase %in% c('PcS','ScP','PcSPcS','ScPScP')){
  if( takeoff<=90){
    maxp=ConvP2Vdepthinv(model$rp-cmbspher,model$vp,model$rp-model$z)
    if( p>maxp){
      
      tt=Inf
      return(list(tt = tt, vdep = vdep, resp = resp))
    }
    vdep=ConvP2Vdepth(p,vp,r,h,model$rp,disconradii) 
    indies=which(vdep<=h) 
    vdep=max(vdep[indies]) 
    if(is.na(model$cmb)){
      
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    }
    if( is.null(vdep)){
      
      vdep=TransformF2Sz(model$cmb+1,model$cmb+1,model$rp)[[1]] 
    }
    p=TransformS2Fp(p,rp)
    dmyz=TransformS2Fz(c(1,1),rp-c(vdep,h),rp)[[2]]
    vdepflat=dmyz[1]
    hflat=dmyz[2]
    if( vdepflat<model$cmb){
      
      tt=Inf
      return(list(tt = tt, vdep = vdep, resp = resp))
    }
    ptt1=CalcTPsum(p,vpflat,zflat,0,model$cmb,0) 
    stt1=CalcTPsum(p,vsflat,zflat,0,model$cmb,0)[[1]] 
    if(phase %in% c('PcS','PcSPcS')){
      tt2=CalcTPsum(p,vpflat,zflat,0,hflat,1)[[1]] 
    }else if(phase %in% c('ScP','ScPScP')){
      tt2=CalcTPsum(p,vsflat,zflat,0,hflat,1) 
    } 
    
    if(phase %in% c('PcS','ScP')){
      tt=ptt1+stt1-tt2
    }else if(phase %in% c('PcSPcS','ScPScP')){
      tt=2*ptt1+2*stt1-tt2
    } 
  }else{
    tt=NaN
  } 
  
  
}else if(phase == 'PKiKP'){
  if( takeoff<=90){
    vdep=ConvP2Vdepth(p,vp,r,h,model$rp,disconradii) 
    indies=which(vdep<=h) 
    vdep=max(vdep[indies]) 
    if( is.na(model$cmb)|is.na(model$icb)){
      
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    } 
    if( is.null(vdep)){
      
      vdep=TransformF2Sz(model$icb+1,model$icb+1,model$rp)[[1]] 
    } 
    p=TransformS2Fp(p,rp)
    dmyz=TransformS2Fz(c(1,1),rp-c(vdep,h),rp)[[2]]
    vdepflat=dmyz[1]
    hflat=dmyz[2]
    if( model$icb-vdepflat>zprecision){
      
      tt=Inf
      return(list(tt = tt, vdep = vdep, resp = resp))
    }
    tt1=CalcTPsum(p,vpflat,zflat,0,model$icb,0)
    tt2=CalcTPsum(p,vpflat,zflat,0,hflat,1) 
    tt=2*tt1-tt2
  }else{
    tt=NaN
  } 
  
  
}else if(phase %in% c('S','SS','SSS')){
  vdep=ConvP2Vdepth(p,vs,r,h,model$rp,disconradii) 
  if( is.na(vdep)&(takeoff<90)){
    
    tt=NaN
    return(list(tt = tt, vdep = vdep, resp = resp))
  }else{
    p=TransformS2Fp(p,rp)
    hflat=TransformS2Fz(c(1,1),rp-h,rp)[[2]]
    if( takeoff<=90){
      indies=which(vdep<h) 
      vdep=max(vdep[indies]) 
      if( vdep<r[length(r)-1]){
        
        tt=Inf
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      if( length(vdep) == 0){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      } 
      if( vdep-(rp-cmbspher)<zprecision){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      
      arriveangle=ConvP2Ang('S',cmbspher,psav,modelsav)
      arriveangle=arriveangle[1] 
      if( isTRUE(arriveangle-90>angleprecision)){
        vdepp=ConvP2Vdepth(p,vp,r,h,model$rp,disconradii)
        indies=which(vdepp<=rp-cmbspher) 
        if( !is.null(indies)){
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
      }
      
    } 
    tt2=CalcTPsum(p,vsflat,zflat,0,hflat,1)
    if( takeoff<=90) {
      vdepflat=TransformS2Fz(c(1,1),rp-vdep,rp)[[2]]
      if( vdepflat-model$cmb>zprecision){
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      tt1=CalcTPsum(p,vsflat,zflat,0,vdepflat,0)
      tt=2*tt1-tt2 
      if(phase =='SS'){
        tt = tt + 2 * tt1
      }else if(phase=='SSS'){
        tt = tt + 4 * tt1
      } 
    }else{
      if(phase=='S'){
        tt=tt2
      }else{
        tt=NaN
      }
    }
  } 
  
}else if(phase %in% c('ScS','ScSScS')){
  if( takeoff<=90){
    
    if( is.na(model$cmb)){
      
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    } 
    
    maxp=ConvP2Vdepthinv(model$rp-cmbspher,model$vs,model$rp-model$z)
    if( p>maxp){
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    } 
    
    
    p=TransformS2Fp(p,rp)
    hflat=TransformS2Fz(c(1,1),rp-h,rp)[[2]]
    
    tt1=CalcTPsum(p,vsflat,zflat,0,model$cmb,0)
    tt2=CalcTPsum(p,vsflat,zflat,0,hflat,1) 
    
    if(phase=='ScS'){
      tt=2*tt1-tt2
    }else if(phase =='ScSScS'){
      tt=4*tt1-tt2
    } 
  }else{
    tt=NaN
  } 
  
}else if(phase %in% c('SKS','SKSSKS','SKKS','SKIKS')){
  if( takeoff<=90){
    vdepsleg=ConvP2Vdepth(p,vs,r,h,model$rp,disconradii)
    indies=which(((vdepsleg<=h)&(vdepsleg<=model$rp-cmbspher))|(is.na(vdepsleg)))
    
    if( is.null(indies)){
      
      tt=Inf
      return(list(tt = tt, vdep = vdep, resp = resp))
    } 
    vdep=ConvP2Vdepth(p,vp,r,model$rp-cmbspher,model$rp,disconradii)
    if( is.na(vdep)){
      
      tt=NaN
      return(list(tt = tt, vdep = vdep, resp = resp))
    }else{
      if(phase %in% c('SKS','SKSSKS','SKKS')){
        if( is.na(model$cmb)){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
        indies=which(vdep<=model$rp-cmbspher) 
      }else if(phase == 'SKIKS'){
        if( is.na(model$cmb)|is.na(model$icb)){
          
          tt=NaN
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
        indies=which(vdep<=model$rp-icbspher) 
      } 
      vdep=max(vdep[indies]) 
      if( vdep<r[length(r)-1]){
        
        tt=Inf
        return(list(tt = tt, vdep = vdep, resp = resp))
      }
      if( is.null(vdep)){
        
        tt=NaN
        return(list(tt = tt, vdep = vdep, resp = resp))
      } 
      p=TransformS2Fp(p,rp)
      dmyz=TransformS2Fz(c(1,1),rp-c(vdep,h),rp)[[2]]
      vdepflat=dmyz[1]
      hflat=dmyz[2]
      
      if(phase %in% c('SKS','SKSSKS','SKKS')){
        if( vdepflat<model$cmb ){
          
          tt=Inf
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
      }else if(phase %in% 'SKIKS'){
        if(vdepflat<=model$icb){
          
          tt=Inf
          return(list(tt = tt, vdep = vdep, resp = resp))
        } 
      } 
      
      tt1=CalcTPsum(p,vsflat,zflat,0,model$cmb,0)
      
      tt2=CalcTPsum(p,vsflat,zflat,0,hflat,1)
      
      ttk=CalcTPsum(p,vpflat,zflat,model$cmb,vdepflat,0)
      
      if(phase %in% c('SKS','SKIKS')){
        tt=2*tt1-tt2+2*ttk
      }else if(phase=='SKSSKS'){
        tt=2*(2*tt1+2*ttk)-tt2
      }else if(phase=='SKKS'){
        tt=2*tt1-tt2+4*ttk
      } 
    } 
  }else{
    tt=NaN
  } 
  
}else{ 
  stop(paste('MKTIM4P: cannot handle phase ', phase, '. Abort.',sep=''))
  tt=NaN
  return(list(tt = tt, vdep = vdep, resp = resp))
} 
  return(list(tt = tt, vdep = vdep, resp = resp))
}

