FindP4Dist <-function(phase,deltalist,h,model,pscan=NULL){
  pres=NULL
  ares=NULL
  dres=NULL
  p=NULL
  a=NULL
  d=NULL
  deltain=NULL
  
  radian=pi/180
  epsilondeg=0.01 
  
  if(is.null(pscan)){
    p=NaN
    pscan=MakePscan(phase,h,model)
  }
  
  if(length(pscan$p)<2){
    return(list(p=NaN,a=NaN,d=NaN,deltain=deltain))
  }
  deltaanz=length(deltalist)
  for(deltacnt in 1:deltaanz){
    delta=deltalist[deltacnt]
    angleanz=length(pscan$angles) 
    for(anglecnt in 1:(angleanz-1)){

      if(pscan$dist[anglecnt]==delta){
        
        
        ares=c(ares, pscan$angles[anglecnt])
        pres=c(pres, pscan$p[anglecnt])
        dres=c(dres, pscan$dist[anglecnt])
        deltain=c(deltain, delta)
      }else{
        if( sign(pscan$dist[anglecnt]-delta)!=sign(pscan$dist[anglecnt+1]-delta) ){
          pad=FindRoots(phase,delta,h,model,c(pscan$angles[anglecnt],pscan$angles[anglecnt+1]),c(pscan$dist[anglecnt],pscan$dist[anglecnt+1]) )
          
          if(!is.na(pad[[3]]) && (pad[[3]] - delta) < 0.1){
            ares=c(ares,pad[[2]])
            pres=c(pres,pad[[1]])
            dres=c(dres,pad[[3]])
            deltain=c(deltain, delta)
          }
          
        }
      }		
    }
  }
  
  if(is.null(pres)){
    p=NaN
  }else{p=pres}
  
  return(list(p=p,a=ares,d=dres,deltain=deltain))
}

