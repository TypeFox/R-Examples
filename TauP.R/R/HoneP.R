HoneP <-function(oldp,oldangle,direction,phase,h,model){
  if(is.na(oldp)){ # break in case of NaN oldp
    return(list(newp=oldp,newangle=oldangle))
  }
  newp=oldp
  newangle=oldangle
  minexponent=-10 
  maxexponent=0
  exponents=minexponent:maxexponent
  crement=10^exponents
  direction=tolower(direction)
  
  if(direction == 'up'){
    crementsign=+1
  }else if(direction == 'down'){
    crementsign=-1
  }else if(direction == 'both'){
    
    aa=HoneP(oldp,oldangle,'up',phase,h,model)
    upp=aa[[1]]
    upangle=aa[[2]]
    
    aa=HoneP(oldp,oldangle,'down',phase,h,model)
    downp=aa[[1]]
    downangle=aa[[2]]
    
    newp=c(upp,downp)
    newangle=c(upangle,downangle)
    nancnt=sum(is.na(newp)) 
    
    if(nancnt == 0){
      deviation=abs(newp-oldp) 
      mindev=min(deviation); 
      indy = which(deviation==mindev)[1]
      newp=newp[indy] 
      newangle=newangle[indy]
    }else if(nancnt == 1){
      
      possible=which(!is.na(newp)) 
      newp=newp[possible]
      newangle=newangle[possible]
    }else if(nancnt == 2){
      
      newp=oldp
      newangle=oldangle
    }else{
      stop('MKFUMBLEP: unexpected number of NaN solutions.')
    }
    
    
    
    return(list(newp=newp,newangle=newangle))
  }else{
    stop(paste('MKFUMBLEP: unknown search direction ', direction))
  } 
  
  done=0
  cnt=1
  maxcnt=length(crement) 
  candidate=oldp
  
  while( done==0 ){
    dist=FindDist4p(phase,h,model,p=candidate)[[1]]
    if( !is.na(dist) ){
      done=1
      newp=candidate
    }else{
      candidate=oldp+crementsign*crement[cnt]
      cnt=cnt+1
      if( cnt==maxcnt){
        done=1
        newp=NaN
      } 
      if( candidate<0){
        done=1
        newp=NaN
      } 
    } 
    
  } 
  
  if( !is.null(newp)){
    
    
    newangle=ConvP2Ang(phase,h,newp,model)
    
    
    if( oldangle<90){
      
      newangle=newangle[length(newangle)]
    }else{
      
      newangle=newangle[1]
    } 
    
  }else{
    
    newangle=NaN
  } 
  
  
  return(list(newp=newp,newangle=newangle))
}

