cosmapfunc=function(cosparamx='CoVol', cosparamy='z', H0=100, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, zrange=c(0,20), step='z', res=100, degen='lo', ref){
  
  paramlistx=c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime','H','OmegaM','OmegaL','OmegaK','Factor','Rate','Sigma8','RhoCrit')
  paramlisty=c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime','H','OmegaM','OmegaL','OmegaK','Factor','Rate','Sigma8','RhoCrit')
  if(! cosparamx %in% paramlistx){stop('cosparamx is not an allowed cosmological parameter, see help options.')}
  if(! cosparamy %in% paramlisty){stop('cosparamy is not an allowed cosmological parameter, see help options.')}
  if(! degen %in% c('lo','hi')){stop('degen option must either be set to lo or hi.')}
  if(cosparamx %in% c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime')){
    pre_x='cosdist'
  }else{
    pre_x='cosgrow'
  }
  if(cosparamy %in% c('z', 'a', 'CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime')){
    pre_y='cosdist'
  }else{
    pre_y='cosgrow'
  }
  
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  
  if(step=='z'){
    zvals=seq(zrange[1],zrange[2],len=res)
  }
  if(step=='logz'){
    zrangelog=log10(1+zrange)
    zvalslog=seq(zrangelog[1],zrangelog[2],len=res)
    zvals=10^zvalslog-1
  }
  if(step=='a'){
    avals=seq(1/(1+zrange[1]),1/(1+zrange[2]),len=res)
    zvals=1/avals-1
  }
  
  if(cosparamx %in% c('CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime','H','RhoCrit')){
    combxparams=list(z=zvals, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
  }
  if(cosparamx %in% c('z', 'a')){
    combxparams=list(z=zvals)
  }
  if(cosparamx %in% c('OmegaM', 'OmegaL', 'OmegaK','Factor')){
    combxparams=list(z=zvals, OmegaM=OmegaM, OmegaL=OmegaL)
  }
  if(cosparamx %in% c('Rate')){
    combxparams=list(z=zvals, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8, fSigma8=fSigma8)
  }
  if(cosparamx %in% c('Sigma8')){
    combxparams=list(z=zvals, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8)
  }
  
  if(cosparamy %in% c('CoDist', 'LumDist', 'AngDist', 'CoDistTran', 'DistMod', 'AngSize', 'CoVol', 'UniAgeAtz','TravelTime','H','RhoCrit')){
    combyparams=list(z=zvals, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL)
  }
  if(cosparamy %in% c('z', 'a')){
    combyparams=list(z=zvals)
  }
  if(cosparamy %in% c('OmegaM', 'OmegaL', 'OmegaK','Factor')){
    combyparams=list(z=zvals, OmegaM=OmegaM, OmegaL=OmegaL)
  }
  if(cosparamy %in% c('Rate')){
    combyparams=list(z=zvals, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8, fSigma8=fSigma8)
  }
  if(cosparamy %in% c('Sigma8')){
    combyparams=list(z=zvals, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8)
  }
  
  tempx=do.call(paste(pre_x,cosparamx,sep=''),combxparams)
  tempy=do.call(paste(pre_y,cosparamy,sep=''),combyparams)
  yatxmin=which.min(tempx)
  yatxmax=which.max(tempx)
  if(degen=='lo'){
    if(yatxmin>1 & yatxmin<res){tempx=tempx[1:yatxmin];tempy=tempy[1:yatxmin]}
    if(yatxmax>1 & yatxmax<res){tempx=tempx[1:yatxmax];tempy=tempy[1:yatxmax]}
  }
  if(degen=='hi'){
    if(yatxmin>1 & yatxmin<res){tempx=tempx[yatxmin:res];tempy=tempy[yatxmin:res]}
    if(yatxmax>1 & yatxmax<res){tempx=tempx[yatxmax:res];tempy=tempy[yatxmax:res]}
  }
  return=approxfun(tempx,tempy)
}

cosmapval=function(val=50, cosparam='CoVol', H0=100, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, zrange=c(-0.99,100), res=100, iter=8, out='cos', degen='lo', ref){
  if(!missing(ref)){
    params=.getcos(ref)
    H0=as.numeric(params['H0'])
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  temp=function(val, cosparam, H0, OmegaM, OmegaL, Sigma8, fSigma8, zlo, zhi, res, iter, out){
    if(cosparam=='DistMod' & zlo==0){zlo=1e-5}
    zrangetemp=c(zlo, zhi)
    for(i in 1:iter){
    tempz=cosmapfunc(cosparamx=cosparam, cosparamy='z', H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8, fSigma8=fSigma8, zrange=zrangetemp, step='z', res=res, degen=degen)
    currentz=tempz(val)
    if(is.na(currentz)){stop('Required cosmological value does not fall within specified redshift/z range')}
    zlonew=max(zrangetemp[1],currentz-(zrangetemp[2]-zrangetemp[1])/res)
    zhinew=currentz+2*(zrangetemp[2]-zrangetemp[1])/res
    zrangetemp=c(zlonew,zhinew)
    }
    if(out=='cos'){
      outdist=unlist(cosdist(currentz, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, age=TRUE, error=T))
      outgrow=unlist(cosgrow(currentz, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8, fSigma8=fSigma8))
      output=c(outdist,outgrow[3:9])
      Error=abs(val-output[[cosparam]])
      if(Error>0){Error=Error/output[[cosparam]]}
      output=c(output,MapError=Error)
    }
    if(out=='z'){
      output=currentz
    }
    return=output
  }
  if(out=='cos'){
    output=as.data.frame(t(Vectorize(temp)(val=val, cosparam=cosparam, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8, fSigma8=fSigma8, zlo=zrange[1], zhi=zrange[2], res=res, iter=iter, out=out)))
  }
  if(out=='z'){
    output=Vectorize(temp)(val=val, cosparam=cosparam, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, Sigma8=Sigma8, fSigma8=fSigma8, zlo=zrange[1], zhi=zrange[2], res=res, iter=iter, out=out)
  }
  return(output)
}