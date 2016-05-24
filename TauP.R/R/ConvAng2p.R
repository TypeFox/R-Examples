ConvAng2p <-
function(phase,h,angle,model=NULL,vp=NULL,vs=NULL,rp=NULL){
radian=pi/180

starttype=tolower(substr(phase,1,1))

if(!is.null(model)){
   
   
   
   focus=InterpModel(model,h,'simple')
   indy=which(focus$z==h)
   vp=focus$vp[indy]
   vs=focus$vs[indy]
   rp=model$rp
   if(!(missing(vp)&missing(vs)&missing(rp))){
   }
}else{
   
   
  if(missing(vp)|missing(vs)|missing(rp)){stop('Must define either model or all of vp, vs, rp')}
}



   if(starttype == 's'){v=vs}
   if(starttype == 'p'){v=vp}   
   if(!(starttype %in% c('p','s'))){stop(paste('MKANGLE2RAYP: phase names cannot start with', substr(phase,1,1)))}

if( length(v)==2 ){
  
  
  whichangle=(angle<=90)+1
}else{
  
  whichangle=angle * 0 + 1}


if( length(v)==2 ){
    vlist=v[whichangle]
}else{
    vlist=v}



angle[is.infinite(angle)]=NaN

p=radian*sin(angle*radian)*(rp-h)/vlist


return(p)}

