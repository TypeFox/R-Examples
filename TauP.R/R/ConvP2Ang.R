ConvP2Ang <-function(phase,h,p,model=NULL,vp=NULL,vs=NULL,rp=NULL){
  radian=pi/180;
  imageps=1e-4;
  starttype=tolower(substr(phase[1],1,1))
  if(!is.null(model)){
    focus=InterpModel(model,h,'simple')
    indy=which(focus$z==h)
    vp=focus$vp[indy]
    vs=focus$vs[indy]
    rp=model$rp
  } 
  
  if(starttype=='s'){
    v=vs
  }else if(starttype=='p'){
    v=vp
  }else{
    stop(paste('ConvP2Ang: phase names cannot start with',starttype))
  } 
  
  if( length(p)>1){
    p=p
  }
  oldwarn = getOption('warn')
  options(warn = -1)
  if(length(v)==1){
    downgoing=asin(0i+p*v/(radian*(rp-h)))/radian
    upgoing=180-downgoing
  }else if(length(v)==2){
    upgoing=180-asin(0i+p*v[1]/(radian*(rp-h)))/radian
    downgoing=asin(0i+p*v[2]/(radian*(rp-h)))/radian
  }else{
    
    stop(paste('ConvP2Ang: overdetermined velocity at z=', h, ' in ', focus$.name, ' (', length(v), 'samples).'))
  } 
  options(warn = oldwarn)
  
  takeoffangle=c(upgoing,downgoing)
  
  indies=which((Re(takeoffangle)==90)&(abs(Im(takeoffangle))<imageps));
  takeoffangle[indies]=Re(takeoffangle[indies])
  indies=which(abs(Im(takeoffangle))>imageps)
  takeoffangle[indies]=NaN
  
  if( h==0){
    takeoffangle=takeoffangle[(length(p)+1):length(takeoffangle)]
  } 
  takeoffangle=Re(takeoffangle)
  return(takeoffangle)
}

