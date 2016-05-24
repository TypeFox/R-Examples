MakePscan <-function(phase,h,imodel){

  ## Pscan: list giving ray parameters and travel distances for a range of takeoff angles.  

### Define internal functions:
  EmptyPscan=function(cnt){
    pscan=list(phase='',
      h=0,
      angles=NULL,
      p=NULL,
      dist=NULL,
      vp=0,
      vs=0,
      starts=NULL,
      ends=NULL)
    return(pscan)
  }
  
  CleanPscan=function(pscan){ ## remove NaNs
    pscannew=pscan
    nandists=which((is.na(pscannew$dist))|(is.infinite(pscannew$dist)))
    if( !is.null(nandists)){
      nonnandists=which((!is.na(pscannew$dist))&(!is.infinite(pscannew$dist)))
      pscannew$dist=pscannew$dist[nonnandists]
      pscannew$angles=pscannew$angles[nonnandists]
      pscannew$p=pscannew$p[nonnandists]
    }
    return(pscannew)
  }
  
  
  UniqueSamples=function(x,y,z){ # filter out repeats
    d=dim(x)
    M = unique(cbind(as.vector(x),as.vector(y),as.vector(z)))
    newx=M[,1]
    newy=M[,2]
    newz=M[,3]
    if(!is.null(d)){
      d[d>1]=length(newx)
      newx=matrix(newx,nrow=d[1],ncol=d[2])
      newy=matrix(newy,nrow=d[1],ncol=d[2])
      newz=matrix(newz,nrow=d[1],ncol=d[2])
    }
    return(list(newx=newx,newy=newy,newz=newz))
  }
  
### Done defining internal functions

  # prepare pscan
  pscan=EmptyPscan(1)
  pscan$phase=phase
  pscan$h=h
  
  angleeps=1e-6 # any differences smaller than this are treated as zero
    
  if(!('criticalrays' %in% names(imodel))){ # if necessary, improve model to specify important ray parameters for discontinuities, etc.
    imodel=ImproveModel(imodel)[[1]]
  }
  
  focus=InterpModel(imodel,h,'simple') # 
  indy=which(focus$z==h)
  pscan$vp=focus$vp[indy]
  pscan$vs=focus$vs[indy]
  
  rayparmlist=NULL
  if('p' %in% tolower(strsplit(phase,'')[[1]]) | 'k' %in% tolower(strsplit(phase,'')[[1]])){ # Have to worry about P waves
    rayparmlist=c(rayparmlist,imodel$criticalrays$p)
  }
  
  if('s' %in% tolower(strsplit(phase,'')[[1]]) | 'j' %in% tolower(strsplit(phase,'')[[1]])){ # Have to worry about S waves
    rayparmlist=c(rayparmlist,imodel$criticalrays$s)
  }
  
  # Have to consider all the takeoff angles that correspond to critical ray parameters
  takeofflist=ConvP2Ang(phase,h,rayparmlist,vp=pscan$vp,vs=pscan$vs,rp=imodel$rp)
  
  if( h!=0 ){ # have to consider both upgoing and downgoing rays--note that both sets of angles are automatically added to takeoff list
    rayparmlist=c(rayparmlist, rayparmlist)
  }
  
  angleinfo=FindPrange(phase,imodel,h,10) 
  
  indies = which((takeofflist > angleinfo$minangle - angleeps) &
    (takeofflist < angleinfo$maxangle + angleeps) &
    (!is.na(takeofflist))) # takeoff angles that are within legal limits for their phase (from FindPrange)
  
  pscan$angles=takeofflist[indies] # update pscan
  pscan$p=rayparmlist[indies]
  
  pscan$angles=c(pscan$angles, angleinfo$minangle, angleinfo$maxangle) # include min/maxangles
  pscan$p=c(pscan$p, ConvAng2p(pscan$phase,pscan$h,c(angleinfo$minangle, angleinfo$maxangle),imodel))
  
  mintemp=HoneP(pscan$p[length(pscan$p)-1],pscan$angles[length(pscan$angles)-1],'both',pscan$phase,pscan$h,imodel) 
  
  maxtemp=HoneP(pscan$p[length(pscan$p)],pscan$angles[length(pscan$angles)],'both',pscan$phase,pscan$h,imodel) 
  
  pscan$p=c(pscan$p, mintemp[[1]], maxtemp[[1]])
  pscan$angles=c(pscan$angles, mintemp[[2]], maxtemp[[2]])
  
  pscan$angles = c(pscan$angles, 90) # add 90 to pscan
  pscan$p = c(pscan$p, ConvAng2p(pscan$phase, pscan$h, 90, imodel))
  
  takeoffanz = length(pscan$angles)
  pscan$dist = FindDist4p(pscan$phase,pscan$h,imodel,takeoff = pscan$angles)[[1]]
    
  sorter=order(pscan$angles, 1:length(pscan$angles), na.last=NA) # sort by angle
  pscan$angles=pscan$angles[sorter]
  pscan$p=pscan$p[sorter]
  pscan$dist=pscan$dist[sorter]
  
  pscan=CleanPscan(pscan) # remove NaNs
  
  intervals=length(pscan$angles)-1
  
  if(intervals < 1){
    return(pscan)
  }
  tt=proc.time()
  
  for( intervalcnt in 1:intervals ){
    
    newpars=OptimizeDist(pscan$angles[intervalcnt+0:1],pscan$dist[intervalcnt+0:1],phase,h,imodel) # look for local maxima/minima in Dist(Angle) over specified interval
    
    
    pscan$angles=c(pscan$angles, newpars[[1]]) 
    pscan$p=c(pscan$p, newpars[[2]]) 
    pscan$dist=c(pscan$dist, newpars[[3]]) 
    
    
  } 
  
  
  pscan=CleanPscan(pscan)
  
  
  sorter=order(pscan$angles)
  pscan$angles = pscan$angles[sorter]
  pscan$p=pscan$p[sorter]
  pscan$dist=pscan$dist[sorter]
  
  aa=UniqueSamples(pscan$angles,pscan$p,pscan$dist)
  pscan$angles=aa[[1]]
  pscan$p=aa[[2]]
  pscan$dist=aa[[3]]
  
  pscan$starts=1:(length(pscan$angles)-1)
  pscan$ends=2:length(pscan$angles)
  
  keep = (pscan$angles - angleinfo[[3]]) < angleeps & (pscan$angles - angleinfo[[2]]) > -angleeps
  pscan$angles = pscan$angles[keep]
  pscan$p = pscan$p[keep]
  pscan$dist = pscan$dist[keep]
  
  return(pscan)
}

