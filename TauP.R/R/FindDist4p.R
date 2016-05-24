FindDist4p <-
function(phase,h,model,p,takeoff)
{
if(missing(p) & missing(takeoff)){stop('mkx4p: must specify one of p, takeoff')}
if(missing(p)){p=ConvAng2p(phase,h,takeoff,model)}
if(missing(takeoff)){takeoff=rep(-Inf,length(p))}
n = max(length(p),length(takeoff))
goodp=rep(TRUE,n)
dist=rep(NaN,n)
segx=as.list(rep(NaN,n))
segz=as.list(rep(NaN,n))
segtyp=as.list(rep(NaN,n))
resp = p
vdep = NULL
radian=pi/180
zprecision=1e-6 
angleprecision=1e-6
imagepsilon=1e-6 
isre = abs(Im(takeoff))<imagepsilon    
goodp = goodp & isre  
takeoff[goodp]=Re(takeoff[goodp])    
goodp = goodp & !is.infinite(p) & !is.na(p) & !is.na(takeoff)
modelsav=model
psav=p
hsav=h
firstchar=substr(phase,1,1)
if((firstchar=='p')|(firstchar=='s')){
   
   remphase=substr(phase,2,nchar(phase)) 
   
   
   if(h==0){
      
      stop('MKX4P: Depth phases do not exist for surface foci!')
    }
      
   dltakeoff=takeoff
   
   dltakeoff[is.infinite(takeoff)]=ConvP2Ang(toupper(firstchar),h,p[is.infinite(takeoff)],model)[1:sum(is.infinite(takeoff))] 
   goodp = goodp & !( (dltakeoff<90) & is.finite(dltakeoff)) & !is.na(dltakeoff)
   a=FindDist4p(toupper(firstchar),h,model,takeoff=dltakeoff)
   dldist=a[[1]]
   dlsegx=a[[2]]
   dlsegz=a[[3]]
   dlsegtyp=a[[4]]
   dlresp=a[[5]]
  
   a=FindDist4p(remphase,0,model,p=p)
   remdist=a[[1]]
   remsegx=a[[2]]
   remsegz=a[[3]]
   remsegtyp=a[[4]]
   remresp=a[[5]]
   bad = is.na(dldist + remdist)
   for(i in which(!bad)){ 
     dist[i]=dldist[i]+remdist[i]
     segx[[i]]=c(dlsegx[[i]],remsegx[[i]][2:length(remsegx[[i]])]+dldist[i])
     segz[[i]]=c(dlsegz[[i]],remsegz[[i]][2:length(remsegz[[i]])])
     segtyp[[i]]=c(dlsegtyp[[i]],remsegtyp[[i]])
     resp[i]=dlresp[i]
   } 
   dist[bad] = NaN
   segx[bad] = NaN
   segz[bad] = NaN
   segtyp[bad] = NaN
   
   return(list(dist=dist,segx=segx, segz=segz, segtyp=segtyp,resp=resp))
 }
a=StripRepetitions(phase)
remphase=a[[1]]
repetitions=a[[2]]
if( repetitions>1 ){
  goodp=goodp & (takeoff<=90 )
  
  basep=rep(NaN,n)
  basep[is.infinite(takeoff)&goodp]=p[is.infinite(takeoff)&goodp]
  
  basep[!is.infinite(takeoff) & goodp]=ConvAng2p(remphase,h,takeoff[!is.infinite(takeoff)],model) 
  
  a=FindDist4p(remphase,h,model,p=basep)
  basedist=a[[1]]
  basesegx=a[[2]]
  basesegz=a[[3]]
  basesegtyp=a[[4]]
  baseresp=a[[5]]
                                        
  
  a=FindDist4p(remphase,0,model,p=basep)
  remdist=a[[1]]
  remsegx=a[[2]]
  remsegz=a[[3]]
  remsegtyp=a[[4]]
  remresp=a[[5]]
  
  dist=basedist+(repetitions-1)*remdist
  for(i in which(goodp)){
    segx[[i]]=basesegx[[i]]
    segz[[i]]=basesegz[[i]]
    segtyp[[i]]=basesegtyp[[i]]
    for(j in 1:(repetitions-1)){
      segx[[i]]=c(segx[[i]],segx[[i]][length(segx[[i]])]+remsegx[[i]][2:length(remsegx[[i]])])
      segz[[i]]=c(segz[[i]],remsegz[[i]][2:length(remsegz[[i]])])
      segtyp[[i]]=c(basesegtyp[[i]],rep(remsegtyp[[i]],(repetitions-1)))
    }
  }
  return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
}
if( max(model$z)<model$rp ){
  anz=length(model$z); 
  model$z=c(model$z,model$rp)
  model$vp=c(model$vp,model$vp[anz])
  model$vs=c(model$vs,model$vs[anz])
  model$rho=c(model$rho,model$rho[anz])
  model$qp=c(model$qp,model$qp[anz])
  model$qs=c(model$qs,model$qs[anz])
}
vp=model$vp
vs=model$vs
rp=model$rp
r=model$rp-model$z
goodp=goodp & p!=0 
goodp = goodp & !is.infinite(p) 
                                        
                                        
h=rp-h
a=TransformS2Fz(vp,rp-r,rp)
vpflat=a[[1]]
zflat=a[[2]]
vsflat=TransformS2Fz(vs,rp-r,rp)[[1]]
cmbspher=model$cmb 
icbspher=model$icb 
alld=c(model$conr, model$moho, model$d410, model$d520, model$d660, model$cmb, model$icb, model$dz)
alld=TransformS2Fz(alld,alld,model$rp)[[2]]
model$conr=alld[1]
model$moho=alld[2]
model$d410=alld[3]
model$d520=alld[4]
model$d660=alld[5]
model$cmb=alld[6]
model$icb=alld[7]
model$dz=alld[8:length(alld)]
if(sum(is.na(model$dz))){model$dz=NULL}
disconradii=FindDiscon(modelsav)
lowermostflat=zflat[length(zflat)-1]
if( is.na(model$cmb) ){
  model$cmb=lowermostflat
  modelsav$cmb=modelsav$z[length(modelsav$z)-1]
  cmbspher=modelsav$z[length(modelsav$z)-1]
}
if( is.na(model$icb) ){
  model$icb=lowermostflat
  modelsav$icb=modelsav$z[length(modelsav$z)-1]
  icbspher=modelsav$z[length(modelsav$z)-1]
}
if(phase == 'PS'){
  aa=FindDist4p('P',hsav,modelsav,p=psav,takeoff=takeoff)
  pdist=aa[[1]]
  psegx=aa[[2]]
  psegz=aa[[3]]
  psegtyp=aa[[4]]
  
  goodp=goodp & (takeoff<=90)
  reflecttakeoff=ConvP2Ang('S',0,psav,modelsav)
  aa=FindDist4p('S',0,modelsav,takeoff=reflecttakeoff)
  sdist=aa[[1]]
  ssegx=aa[[2]]
  ssegz=aa[[3]]
  ssegtyp=aa[[4]]
  goodp=goodp & !is.na(sdist + pdist)
  for( i in which(goodp) ){
    if( ((max(ssegz[[i]])-modelsav$cmb)<zprecision)&((max(psegz[[i]])-modelsav$cmb)<zprecision) ){
      
      dist[i]=pdist[i]+sdist[i]
      segx[[i]]=c(psegx[[i]], psegx[[i]][length(psegx[[i]])]+ssegx[[i]][2:length(ssegx[[i]])])	
      segz[[i]]=c(psegz[[i]], ssegz[[i]][2:length(ssegz[[i]])])
      if( (!is.na(pdist[i]))&(!is.na(sdist[i]))&(is.finite(pdist[i]))&(is.finite(sdist[i])) ){
        segtyp[[i]]=c(psegtyp[[i]],ssegtyp[[i]])}
    }
  } 
  
  
}else if(phase == 'SP'){
  aa=FindDist4p('S',hsav,modelsav,p=psav,takeoff=takeoff)
  sdist=aa[[1]]
  ssegx=aa[[2]]
  ssegz=aa[[3]]
  ssegtyp=aa[[4]]
  goodp = goodp & ( takeoff<=90 )
  reflecttakeoff=ConvP2Ang('P',0,psav,modelsav)
  aa=FindDist4p('P',0,modelsav,p=psav,takeoff=reflecttakeoff)
  pdist=aa[[1]]
  psegx=aa[[2]]
  psegz=aa[[3]]
  psegtyp=aa[[4]]
  goodp = goodp & !is.na(pdist) & !is.na(sdist)
  for( i in which(goodp)){
    if( (abs(max(psegz[[i]])-modelsav$cmb)>zprecision)&(abs(max(ssegz[[i]])-modelsav$cmb)>zprecision) ){
      dist[i]=pdist[i]+sdist[i]
      segx[[i]]=c(ssegx[[i]],ssegx[[i]][length(ssegx[[i]])]+psegx[[i]][2:length(psegx[[i]])])
      segz[[i]]=c(ssegz[[i]],psegz[[i]][2:length(psegz[[i]])])
      if( (!is.na(pdist[i]))&(!is.na(sdist[i]))&(!is.infinite(pdist[i]))&(!is.infinite(sdist[i])) ){
        segtyp[[i]]=c(ssegtyp[[i]],psegtyp[[i]])
      }
    }
  }
  
  
      } else if(phase %in% c('P','PP','PPP','PKP','PKPPKP','PKIKP','PKIKPPKIKP')){
        
                                        
                                        
                                        
                                        
                                        
        time=proc.time()
        vdep=rep(NaN,n)
        for(i in which(goodp)){
          vdep[i]=ConvP2Vdepth(p[i],vp,r,h,model$rp,disconradii) 
        }
                                        
                                        
        goodp = goodp & !(is.na(vdep) & takeoff<90)
        p=TransformS2Fp(p,rp)
        hflat=TransformS2Fz(c(1,1),rp-h,rp)[[2]]
        indies=which(vdep<=h) 
                                        
                                        
                                        
        dist[ vdep<r[length(r)-1] ]=Inf
        goodp = goodp & !(vdep<r[length(r)-1])
                                        
                                        
        dmyz=TransformS2Fz(rep(1,length(vdep)+1),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[2:length(dmyz)]
        hflat=dmyz[1]
        
        if( phase %in% c('PKP','PKPPKP') ){
          if( is.na(model$cmb) ){
            print('Warning: CMB undefined!')	
            return(list(dist=p+NaN,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
          }
          goodp = goodp & !((model$cmb-vdepflat>zprecision)|(model$icb-vdepflat<zprecision) )
          
        }
        if( phase %in% c('PKIKP','PKIKPPKIKP') ){
          if( is.na(model$cmb)|is.na(model$icb) ){
            print('Warning: CMB or ICB is undefined')
            return(list(dist=NaN,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
          }
          goodp = goodp & !( model$icb-vdepflat>zprecision )
        }
        if( phase == 'P' ){
          goodp = goodp & !(vdepflat-model$cmb>zprecision & takeoff<90) 
        }
          
        if( phase %in% c('PP','PPP') ){
          goodp = goodp & !(vdepflat-model$cmb>zprecision)
          
        }
        
                                        
        dist2=NULL
        segx2=list()
        segz2=list()
        for(i in which(goodp)){   
          aa=CalcXPsum(p[i],vpflat,zflat,0,hflat,1) 
          dist2[i]=aa[[1]]
          segx2[[i]]=aa[[2]]
          segz2[[i]]=aa[[3]]
        }
                                        
        dist1=NULL
        segx1=list()
        segz1=list()
        segz1s=list()
        segz2s=list()
        segxs=list()
        segx2s=list()
        for(i in which(goodp & takeoff<=90) ){
          
          aa=CalcXPsum(p[i],vpflat,zflat,0,vdepflat[i],0)
          dist1=aa[[1]]
          segx1=aa[[2]]
          segz1=aa[[3]]
                                        
          dist[i]=TransformF2Sdist(2*dist1-dist2[i],rp)
          
          segz1=TransformF2Sz(segz1,segz1,rp)[[2]] 
          segz2s=TransformF2Sz(segz2[[i]],segz2[[i]],rp)[[2]] 
          segxs=cumsum(TransformF2Sdist(segx1,rp))  
          segx2s=cumsum(TransformF2Sdist(segx2[[i]],rp)); 
          
          indies=length(segxs):1
          segx[[i]]=c(segxs,2*segxs[length(segxs)]-segxs[indies])
          segz[[i]]=c(segz1,segz1[indies])
          
          
          if( h!=rp ){
            indies=c(which(segz[[i]][1:length(segz1)]>(rp-h)),(length(segz1)+1):length(segz[[i]]) ); 
            segx[[i]]=segx[[i]][indies]
            segz[[i]]=segz[[i]][indies]
            
            segx[[i]]=c(0,segx[[i]]-segx2s[length(segx2s)])
            segz[[i]]=c(rp-h, segz[[i]])
          }
          
          
          if( phase %in% c('PP','PKPPKP','PKIKPPKIKP') ){
            indies=length(segxs):1
            segxss=c(segxs,2*max(segxs)-segxs[indies])
            segzss=c(segz1,segz1[indies])
            segx[[i]]=c(segx[[i]],segxss+dist[i])
            segz[[i]]=c(segz[[i]],segzss)
            dist[i]=dist[i]+2*TransformF2Sdist(dist1,rp)
          }
          if( phase == 'PPP' ){
            indies=length(segxs):1
            dist1=TransformF2Sdist(dist1,rp)
            segxss=c(segxs, 2*max(segxs)-segxs[indies])+dist[i]
            segzss=c(segz1, segz1[indies])
            segx[[i]]=c(segx[[i]], segxss, segxss+2*dist1)
            segz[[i]]=c(segz[[i]], segzss, segzss)
            dist[i]=dist[i] + 4*dist1
          }
          segtyp[[i]]=rep('P',length(segx[[i]])-1)
        }
        if( phase == 'P' ){
          for(i in which(goodp & (takeoff>90))){
            
            segz2s=TransformF2Sz(segz2[[i]],segz2[[i]],rp)[[2]]
            segx2s=cumsum(TransformF2Sdist(segx2[[i]],rp)); 
            dist[i]=TransformF2Sdist(dist2[[i]],rp)
            segx[[i]]=max(segx2s)-rev(segx2s)
            segz[[i]]=rev(segz2s)
            segtyp[[i]]=rep('P',length(segx[[i]])-1)
          }
        }
                                        
        
        
      } else if(phase == 'PKKP'){
        goodp = goodp & (takeoff<=90 )
        for(i in which(goodp)){
          vdep[i]=ConvP2Vdepth(p[i],vp,r,h,model$rp,disconradii); 
        }
        goodp = goodp & (vdep<=h)
        
        
        if( is.na(model$cmb) ){
          return(list(dist=NaN,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
        }
        goodp = goodp & !( vdep<r[length(r)-1] )
        
        p=TransformS2Fp(p,rp)
        dmyz=TransformS2Fz(rep(1,length(vdep)+1),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[2:length(dmyz)]
        hflat=dmyz[1]
        
        goodp = goodp & !(vdepflat<=model$cmb) 
        goodp = goodp & !(vdepflat>=model$icb) 
        
        for(i in which(goodp)){
            aa=CalcXPsum(p[i],vpflat,zflat,0,model$cmb,0)
            dist1=aa[[1]]
            segx1=aa[[2]]
            segz1=aa[[3]]
            aa=CalcXPsum(p[i],vpflat,zflat,0,hflat,1); 
            dist2=aa[[1]]
            segx2=aa[[2]]
            segz2=aa[[3]]
            aa=CalcXPsum(p[i],vpflat,zflat,model$cmb,vdepflat[i],0);
            kdist=aa[[1]]
            ksegx=aa[[2]]
            ksegz=aa[[3]]
            
          
            segz1s=TransformF2Sz(segz1,segz1,rp)[[2]]
            segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
            ksegzs=TransformF2Sz(ksegz,ksegz,rp)[[2]] 
            segxs=cumsum(TransformF2Sdist(segx1,rp))  
            segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
            ksegxs=cumsum(TransformF2Sdist(ksegx,rp)) 
            
            indies=length(segxs):1
            kindies=length(ksegxs):1
            klegx=c(ksegxs, 2*max(ksegxs)-ksegxs[kindies])
            klegz=c(ksegzs,ksegzs[kindies])
            
            segx[[i]]=c(segxs, klegx+segxs[length(segxs)],klegx+segxs[length(segxs)]+klegx[length(klegx)],2*segxs[length(segxs)]-segxs[indies]+2*klegx[length(klegx)])
            segz[[i]]=c(segz1s,klegz,klegz,segz1s[indies])
            dist[i]=TransformF2Sdist(2*dist1-dist2+4*kdist,rp)
            
            
            if( h!=rp ){
              indies=c(which(segz[[i]][1:length(segz1s)]>(rp-h)),(length(segz1s)+1):length(segz[[i]])) 
              segx[[i]]=segx[[i]][indies]
              segz[[i]]=segz[[i]][indies]
              
              segx[[i]]=c(0,segx[[i]]-segx2s[length(segx2s)])
              segz[[i]]=c(rp-h,segz[[i]])
          }
          segtyp[[i]]=rep('P',length(segx[[i]])-1)
          indies=(segz[[i]]>=cmbspher) 
          segtyp[[i]][indies]='P' 
          indies=(segtyp[[i]][1:(length(segtyp)-1)]!=segtyp[2:length(segtyp)]) 
          segtyp[[i]][indies]='P' 
          } 
        
      } else if( phase %in% c('PcP','PcPPcP') ){ 
        goodp = goodp & (takeoff <= 90)
        maxp=ConvP2Vdepthinv(model$rp-cmbspher,model$vp,model$rp-model$z)
        dist[goodp & p>maxp]=Inf
        goodp = goodp & p<=maxp
        vdep=rep(NaN,length(p))
        for(i in which(goodp)){
          vdep[i]=ConvP2Vdepth(p[i],vp,r,h,model$rp,disconradii) 
        }
        goodp = goodp & (vdep<=h)
        vdep[!(vdep<=h)]=NaN
        if( is.na(model$cmb) ){
          dist=NaN
          warning('model$cmb = NaN')
          return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
        } 
          
        p=TransformS2Fp(p,rp)
        dmyz=TransformS2Fz(rep(1,length(vdep)+1),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[2:length(dmyz)]
        hflat=dmyz[1]
        dist[vdepflat<model$cmb]=Inf 
        goodp = goodp & vdepflat>model$cmb
          
        for(i in which(goodp)){
          a=CalcXPsum(p[i],vpflat,zflat,0,model$cmb,0)
          dist1=a[[1]]
          segx1=a[[2]]	
          segz1=a[[3]]
          a=CalcXPsum(p[i],vpflat,zflat,0,hflat,1)
          dist2=a[[1]]
          segx2=a[[2]]
          segz2=a[[3]]
          
          
          segz1s=TransformF2Sz(segz1,segz1,rp)[[2]]
          segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
          segxs=cumsum(TransformF2Sdist(segx1,rp))  
          segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
          
          
          if( phase == 'PcP' ){
            dist[i]=TransformF2Sdist(2*dist1-dist2,rp)
            indies=length(segxs):1
            segx[[i]]=c(segxs, 2*segxs[length(segxs)]-segxs[indies])
            segz[[i]]=c(segz1s, segz1s[indies])
          }
          if( phase == 'PcPPcP' ){
            dist[i]=TransformF2Sdist(4*dist1-dist2,rp)
            indies=length(segxs):1
            segx[[i]]=c(segxs, 2*segxs[length(segxs)]-segxs[indies])	
            segx[[i]]=c(segx[[i]], segx[[i]]+segx[[i]][length(segx[[i]])])
            segz[[i]]=c(segz1s, segz1s[indies], segz1s, segz1s[indies])
          }
          
          if( h!=rp ){
            indies=c(which(segz[[i]][1:length(segz1s)]>(rp-h)),(length(segz1s)+1):length(segz[[i]])); 
            segx[[i]]=segx[[i]][indies]
            segz[[i]]=segz[[i]][indies]
            
            segx[[i]]=c(0, segx[[i]]-segx2s[length(segx2s)])
            segz[[i]]=c(rp-h, segz[[i]])
          }
          segtyp[[i]]=rep('P',length(segx[[i]]))
        }
        
      } else if( phase %in% c('PcS','ScP','PcSPcS','ScPScP') ){ 
        goodp = goodp & takeoff<=90
        maxp=ConvP2Vdepthinv(model$rp-cmbspher,model$vp,model$rp-model$z)
        dist[p>maxp]=Inf
        goodp = goodp & p<=maxp
        vdep = rep(NaN, length(p))
        for(i in which(goodp)){
          vdep[i]=ConvP2Vdepth(p[i],vp,r,h,model$rp,disconradii) 
        }
        if( is.na(model$cmb) ){
          dist=NaN
          warning('model$cmb=NaN')
          return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
        } 
        p=TransformS2Fp(p,rp)
        dmyz=TransformS2Fz(rep(1,length(vdep)+1),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[2:length(dmyz)]
        hflat=dmyz[1]
        dist[vdepflat<model$cmb]=Inf 
        goodp = goodp & (vdepflat>=hflat) & (vdepflat >= model$cmb)
        
        for(i in which(goodp)){
          aa=CalcXPsum(p[i],vpflat,zflat,0,model$cmb,0) 
          pdist1=aa[[1]]
                                        
          psegx1=aa[[2]]
          psegz1=aa[[3]]
          aa=CalcXPsum(p[i],vsflat,zflat,0,model$cmb,0); 
          sdist1=aa[[1]]
          ssegx1=aa[[2]]
          ssegz1=aa[[3]]
          if( phase %in% c('PcS','PcSPcS') ){
            aa=CalcXPsum(p[i],vpflat,zflat,0,hflat,1); 
            dist2=aa[[1]]
            segx2=aa[[2]]
            segz2=aa[[3]]
          }  
          if( phase %in% c('ScP','ScPScP') ){
            aa=CalcXPsum(p[i],vsflat,zflat,0,hflat,1); 
            dist2=aa[[1]]
            segx2=aa[[2]]
            segz2=aa[[3]]
          } 
          
          
          if( phase %in% c('PcS','ScP') ){
            dist[i]=TransformF2Sdist(pdist1+sdist1-dist2,rp)
          }   
          if( phase %in% c('PcSPcS','ScPScP') ){
            dist[i]=TransformF2Sdist(2*pdist1+2*sdist1-dist2,rp)
          } 
          
          psegz1s=TransformF2Sz(psegz1,psegz1,rp)[[2]] 
          segz2=TransformF2Sz(segz2,segz2,rp)[[2]] 
          psegxs=cumsum(TransformF2Sdist(psegx1,rp))  
          ssegz1s=TransformF2Sz(ssegz1,ssegz1,rp)[[2]] 
          psegxs=cumsum(TransformF2Sdist(psegx1,rp))  
          ssegxs=cumsum(TransformF2Sdist(ssegx1,rp))  
          segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
          
          indies=(length(psegxs):1);  
          if( phase == 'PcS'){
            segx[[i]]=c(psegxs,psegxs[length(psegxs)]+ssegxs[length(ssegxs)]-ssegxs[indies[2:length(indies)]])
            segz[[i]]=c(psegz1s, ssegz1s[indies[2:length(indies)]])
            segtyp[[i]]=c(rep('P',length(psegxs)-1), rep('S',length(indies)-1) )
          }
          if( phase == 'PcSPcS'){
            segx[[i]]=c(psegxs, psegxs[length(psegxs)]+ssegxs[length(ssegxs)]-ssegxs[indies])
            segx[[i]]=c(segx[[i]],segx[[i]]+segx[[i]][length(segx[[i]])])
            segz[[i]]=c(psegz1s, ssegz1s[indies], psegz1s, ssegz1s[indies])
            segtyp[[i]]=c(rep('P',length(psegxs)),
              rep('S',length(ssegxs)),
              rep('P',length(psegxs)),
              rep('S',length(ssegxs)-1) )
          }	  
          if( phase == 'ScP'){
            segx[[i]]=c(ssegxs,ssegxs[length(ssegxs)]+psegxs[length(psegxs)]-psegxs[indies[2:length(indies)]])
            segz[[i]]=c(ssegz1s, psegz1s[indies[2:length(indies)]])
            segtyp[[i]]=c(rep('S',length(ssegxs)-1), rep('P',length(indies)-1) )
		}   
          if( phase == 'ScPScP' ){
            segx[[i]]=c(ssegxs, ssegxs[length(ssegxs)]+psegxs[length(psegxs)]-psegxs[indies])
            segx[[i]]=c(segx[[i]],segx[[i]]+segx[[i]][length(segx[[i]])])
            segz[[i]]=c(ssegz1s,psegz1s[indies],ssegz1s,psegz1s[indies])
            segtyp[[i]]=c(rep('S',length(ssegxs)),
              rep('P',length(psegxs)),
              rep('S',length(ssegxs)),
              rep('P',length(psegxs)-1) )
          } 
          
          if( h!=rp ){
            segzlen=length(segz[[i]])
            if( phase %in% c('PcS','PcSPcS') ){
              indies=c(which(segz[[i]][1:length(psegz1s)]>(rp-h)),(length(psegz1s)+1):length(segz[[i]]) ) 
              removed=segzlen-length(indies) 
            }
            if( phase %in% c('ScP','ScPScP') ){
              indies=c(which(segz[[i]][1:length(ssegz1s)]>(rp-h)),(length(ssegz1s)+1):length(segz[[i]]) ); 
              removed=segzlen-length(indies) 
            }
            segx[[i]]=segx[[i]][indies]
            segz[[i]]=segz[[i]][indies]
            segtyp[[i]]=segtyp[[i]][removed:length(segtyp[[i]])]
            
            segx[[i]]=c(0,segx[[i]]-segx2s[length(segx2s)])
            segz[[i]]=c(rp-h,segz[[i]])
          } 
        }
      } else if( phase == 'PKiKP' ){ 
        goodp = goodp & ( takeoff<=90 )
          for(i in which(goodp)){
            vdep[i]=ConvP2Vdepth(p[i],vp,r,h,model$rp,disconradii) 
          }
          goodp = goodp & (vdep<=h) 
             if( is.na(model$icb)|is.na(model$cmb) ){
               dist=NaN+p 
               return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
             } 
        vdep[is.na(vdep)]=TransformF2Sz(model$icb+1,model$icb+1,model$rp)[[2]] 
        p=TransformS2Fp(p,rp)
        dmyz=TransformS2Fz(rep(1,length(vdep)+1),rp-c(h,vdep),rp)
        vdepflat=dmyz[[2]][1:length(vdep)+1]
        hflat=dmyz[[2]][1]
        dist[model$icb-vdepflat>zprecision]=Inf
        goodp = goodp & !(model$icb-vdepflat>zprecision)
        for( i in which(goodp) ){
          aa=CalcXPsum(p[i],vpflat,zflat,0,model$icb,0)
          dist1=aa[[1]]
          segx1=aa[[2]]
          segz1=aa[[3]]
          aa=CalcXPsum(p[i],vpflat,zflat,0,hflat,1) 
          dist2=aa[[1]]
          segx2=aa[[2]]
          segz2=aa[[3]]
          dist[i]=TransformF2Sdist(2*dist1-dist2,rp)
          
          segz1s=TransformF2Sz(segz1,segz1,rp)[[2]]
          segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
          segxs=cumsum(TransformF2Sdist(segx1,rp))  
          segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
          
          indies=(length(segxs):1)
          segx[[i]]=c(segxs,2*segxs[length(segxs)]-segxs[indies])
          segz[[i]]=c(segz1s,segz1s[indies])
          
          if( h!=rp ){
            indies=c(which(segz[[i]][1:length(segz1s)]>(rp-h)),(length(segz1s)+1):length(segz[[i]])) 
            segx[[i]]=segx[[i]][indies]
            segz[[i]]=segz[[i]][indies]
            
            segx[[i]]=c(0,segx[[i]]-segx2s[length(segx2s)])
            segz[[i]]=c(rp-h, segz[[i]])
          } 
          segtyp[[i]]=rep('P',length(segx[[i]])-1)
        } 
        
      } else if( phase %in% c('S','SS','SSS') ){ 
        for(i in which(goodp)){
          vdep[i]=ConvP2Vdepth(p[i],vs,r,h,model$rp,disconradii); 
        }
        vdep[takeoff==90] = h
                                        
                                        
                                        
        
        dist[is.na(vdep)&(takeoff<90)]=NaN
        goodp = goodp & !(is.na(vdep)&(takeoff<90))
        p=TransformS2Fp(p,rp)
        hflat=TransformS2Fz(1,rp-h,rp)[[2]]
        vdep[!(vdep <= h)] = NaN
        dist[is.na(vdep)]=NaN
        dist[ vdep<r[length(r)-1] ] = Inf
        dist[ vdep<rp-cmbspher ] = NaN
            
        dmyz=TransformS2Fz(rep(1,length(vdep)+1),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[2:length(dmyz)]
        hflat=dmyz[1]
        if(phase == 'S'){
          goodp = goodp & !( vdepflat>model$cmb & takeoff<90) 
        } else {
          goodp = goodp & !( vdepflat>model$cmb ) 
        }
        for(i in which(goodp)){
          aa=CalcXPsum(p[i],vsflat,zflat,0,hflat,1)
          dist2=aa[[1]]
          segx2=aa[[2]]
          segz2=aa[[3]]
          if( takeoff[i]<=90 ){ 
            aa=CalcXPsum(p[i],vsflat,zflat,0,vdepflat[i],0)
	    dist1=aa[[1]]
	    segx1=aa[[2]]
	    segz1=aa[[3]]
            dist[i]=TransformF2Sdist(2*dist1-dist2,rp)
            
            segz1s=TransformF2Sz(segz1,segz1,rp)[[2]]
            segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
            segxs=cumsum(TransformF2Sdist(segx1,rp))  
            segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
            
            indies=(length(segxs):1)
            segx[[i]]=c(segxs,2*segxs[length(segxs)]-segxs[indies])
            segz[[i]]=c(segz1s, segz1s[indies])
            
            if( h!=rp ){
              indies=c(which(segz[[i]][1:length(segz1s)]>(rp-h)), (length(segz1s)+1):length(segz[[i]])) 
              segx[[i]]=segx[[i]][indies]
              segz[[i]]=segz[[i]][indies]
              
              segx[[i]]=c(0, segx[[i]]-segx2s[length(segx2s)])
              segz[[i]]=c(rp-h,segz[[i]])
            } 
            if( phase == 'SS' ){
              indies=(length(segxs):1)
              segxss=c(segxs,2*segxs[length(segxs)]-segxs[indies])
              segzss=c(segz1s,segz1s[indies])
              segx[[i]]=c(segx[[i]],segxss+dist[i])
              segz[[i]]=c(segz[[i]],segzss)
              dist[i]=dist[i]+2*TransformF2Sdist(dist1,rp)
            }  
            if( phase == 'SSS' ){
              indies=(length(segxs):1)
              dist1=TransformF2Sdist(dist1,rp)
              segxss=c(segxs,2*segxs[length(segxs)]-segxs[indies])+dist[i]
              segzss=c(segz1s,segz1s[indies])
              segx[[i]]=c(segx[[i]],segxss,segxss+2*dist1)
              segz[[i]]=c(segz[[i]],segzss,segzss)
              dist[i]=dist[i]+4*dist1
            }
            segtyp[[i]]=rep('S',length(segx[[i]])-1)
          }else{
            if( phase == 'S'){
              segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
              segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
              dist[i]=TransformF2Sdist(dist2,rp)
              segx[[i]]=segx2s[length(segx2s)]-rev(segx2s)
              segz[[i]]=rev(segz2s)
              segtyp[[i]]=rep('S',length(segx[[i]])-1)
            }else{
              dist[i]=NaN
            } 
          } 
        } 
        
      }else if( phase %in% c('ScS','ScSScS') ){ 
        goodp = goodp & takeoff<=90
        for(i in which(goodp)){
          vdep[i]=ConvP2Vdepth(p[i],vs,r,h,model$rp,disconradii) 
        }
        dist[!is.null(vdep) && !is.na(vdep) & vdep>(rp-cmbspher+zprecision)] = Inf
        if( is.na(model$cmb) ){
          dist=NaN + p
          return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
        }    
        if(!is.null(vdep)){
          vdep[is.na(vdep)]=TransformF2Sz(+1,+1,model$rp)[[2]] 
        }
        
        p=TransformS2Fp(p,rp)
        dmyz=TransformS2Fz(rep(1,1+length(p)),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[1:length(p)+1]
        hflat=dmyz[1]
        for(i in which(goodp)){
          aa=CalcXPsum(p[i],vsflat,zflat,0,model$cmb,0)
          dist1=aa[[1]]
          segx1=aa[[2]]
          segz1=aa[[3]]
          aa=CalcXPsum(p[i],vsflat,zflat,0,hflat,1); 
          dist2=aa[[1]]
          segx2=aa[[2]]	
          segz2=aa[[3]]
          if( phase == 'ScS' ){
            dist[i]=TransformF2Sdist(2*dist1-dist2,rp)
          }
          if( phase == 'ScSScS' ){
            dist[i]=TransformF2Sdist(4*dist1-dist2,rp)
          }
          
          segz1s=TransformF2Sz(segz1,segz1,rp)[[2]]
          segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
          segxs=cumsum(TransformF2Sdist(segx1,rp));  
          segx2s=cumsum(TransformF2Sdist(segx2,rp)); 
          
          if( phase == 'ScS' ){
            indies=(length(segxs):1)
            segx[[i]]=c(segxs,2*segxs[length(segxs)]-segxs[indies])
            segz[[i]]=c(segz1s,segz1s[indies])
          }   
          if( phase == 'ScSScS' ){
            indies=(length(segxs):1);
            segx[[i]]=c(segxs,2*segxs[length(segxs)]-segxs[indies])
            segx[[i]]=c(segx[[i]],segx[[i]]+segx[[i]][length(segx[[i]])])
            segz[[i]]=c(segz1s,segz1s[indies],segz1s,segz1s[indies])
          }
          
          if( h!=rp ){
            indies=c(which(segz[[i]][1:length(segz1s)]>(rp-h)),(length(segz1s)+1):length(segz[[i]])) 
            segx[[i]]=segx[[i]][indies]
            segz[[i]]=segz[[i]][indies]
            
            segx[[i]]=c(0,segx[[i]]-segx2s[length(segx2s)])
            segz[[i]]=c(rp-h,segz[[i]])
          } 
          segtyp[[i]]=rep('S',length(segx[[i]])-1)
        }
      }else if( phase %in% c('SKS','SKSSKS','SKKS','SKIKS') ){
        goodp = goodp & takeoff<=90
        vdep=rep(NaN,length(p))
        vdepsleg=rep(NaN,length(vdep))
        for(i in which(goodp)){
          vdepsleg[i]=ConvP2Vdepth(p[i],vs,r,h,model$rp,disconradii)
          vdep[i]=ConvP2Vdepth(p[i],vp,r,model$rp-cmbspher,model$rp,disconradii)
        }
        goodp = goodp & (((vdepsleg<=h)&(vdepsleg<=model$rp-cmbspher))|(is.na(vdepsleg)))
        dist[!(((vdepsleg<=h)&(vdepsleg<=model$rp-cmbspher))|(is.na(vdepsleg)))]=Inf
        if( phase %in% c('SKS','SKSSKS','SKKS') ){
          if( is.na(model$cmb) ){
            dist=NaN
            return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
          } 
          goodp=goodp & vdep <= (model$rp-cmbspher) 
          goodp=goodp & vdep > (model$rp-icbspher) 
        }   
        if( phase == 'SKIKS' ){
          if( is.na(model$cmb)|is.na(model$icb) ){
            dist=NaN
            return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
          } 
          goodp=goodp & vdep<=model$rp-icbspher 
        }
        dist[is.na(vdep)]=NaN
        goodp = goodp & !is.na(vdep)
        goodp = goodp & !( vdep < r[length(r)-1] )
        dist[vdep < r[length(r)-1] ] = Inf
        p=TransformS2Fp(p,rp)
        dmyz=TransformS2Fz(rep(1,1+length(vdep)),rp-c(h,vdep),rp)[[2]]
        vdepflat=dmyz[2:length(dmyz)]
        hflat=dmyz[1]
        if( phase %in% c('SKS','SKSSKS','SKKS') ){
          goodp = goodp & !(vdepflat < model$cmb )
          dist[vdepflat < model$cmb]=NaN
        }     
        if( phase == 'SKIKS' ){
          if( is.na(model$icb) ){
            dist=NaN
            return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
          } 
          goodp = goodp & !( vdepflat<=model$icb )
          dist[vdepflat<=model$icb] = NaN
        }
        for(i in which(goodp)){
          aa=CalcXPsum(p[i],vsflat,zflat,0,model$cmb,0)
          dist1=aa[[1]]
          segx1=aa[[2]]
          segz1=aa[[3]]
          aa=CalcXPsum(p[i],vsflat,zflat,0,hflat,1); 
          dist2=aa[[1]]
          segx2=aa[[2]]
          segz2=aa[[3]]
          aa=CalcXPsum(p[i],vpflat,zflat,model$cmb,vdepflat[i],0);
          kdist=aa[[1]]
          ksegx=aa[[2]]
          ksegz=aa[[3]]
          
          segz1s=TransformF2Sz(segz1,segz1,rp)[[2]]
          segz2s=TransformF2Sz(segz2,segz2,rp)[[2]]
          ksegzs=TransformF2Sz(ksegz,ksegz,rp)[[2]] 
          segxs=cumsum(TransformF2Sdist(segx1,rp))  
          segx2s=cumsum(TransformF2Sdist(segx2,rp)) 
          ksegxs=cumsum(TransformF2Sdist(ksegx,rp)) 
          
          indies=(length(segxs):1)
          kindies=(length(ksegxs):1)
          klegx=c(ksegxs,2*ksegxs[length(ksegxs)]-ksegxs[kindies])
          klegz=c(ksegzs,ksegzs[kindies])
          
          if( phase %in% c('SKS','SKIKS') ){
            segx[[i]]=c(segxs,klegx+segxs[length(segxs)],2*segxs[length(segxs)]-segxs[indies]+klegx[length(klegx)])
            segz[[i]]=c(segz1s,klegz,segz1s[indies])
            dist[i]=TransformF2Sdist(2*dist1-dist2+2*kdist,rp)
          }
          if( phase == 'SKKS' ){
            segx[[i]]=c(segxs,klegx+segxs[length(segxs)],klegx+segxs[length(segxs)]+klegx[length(klegx)],2*segxs[length(segxs)]-segxs[indies]+2*klegx[length(klegx)])
            segz[[i]]=c(segz1s,klegz,klegz,segz1s[indies])
            dist[i]=TransformF2Sdist(2*dist1-dist2+4*kdist,rp);
          }
          if( phase == 'SKSSKS' ){
            segx[[i]]=c(segxs,klegx+segxs[length(segxs)],2*segxs[length(segxs)]-segxs[indies]+klegx[length(klegx)])
            segx[[i]]=c(segx[[i]],segx[[i]]+segx[[i]][length(segx[[i]])])
            segz[[i]]=c(segz1s,klegz,segz1s[indies])
            segz[[i]]=c(segz[[i]],segz[[i]])
            dist[i]=TransformF2Sdist(4*dist1-dist2+4*kdist,rp)
          } 
          
          
          if( h!=rp ){
            indies=c(which(segz[[i]][1:length(segz1s)]>(rp-h)),(length(segz1s)+1):length(segz[[i]])) 
            segx[[i]]=segx[[i]][indies]
            segz[[i]]=segz[[i]][indies]
            
            segx[[i]]=c(0,segx[[i]]-segx2s[length(segx2s)])
            segz[[i]]=c(rp-h,segz[[i]])
          } 
          segtyp[[i]]=rep('S',length(segx[[i]])-1)
          indies=which(segz[[i]]>cmbspher) 
          segtyp[[i]][indies]='P' 
          indies=which(segtyp[[i]][2:length(segtyp[[i]])-1] != segtyp[[i]][2:length(segtyp[[i]])]) 
          segtyp[[i]][indies]='S' 
        }
      } else { 
        stop(paste('MKX4P: cannot handle phase ', phase, '. Returning NaN.',sep=''))
        dist=NaN
        return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
      } 
bad = is.na(dist) | is.infinite(dist)
segx[bad]=NaN
segz[bad]=NaN
segtyp[bad]=NaN
return(list(dist=dist,segx=segx,segz=segz,segtyp=segtyp,resp=resp))
}

