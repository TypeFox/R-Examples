ReadND <-
function(filename,verbose=TRUE){

# turn off warnings--coercion creates NAs, which is normal here
  options(warn=-1)
  
  vmod=-1     
  dz=NULL     
  dname=NULL  
  model=-1    
  
  pradius=NULL
  pubyear=NULL
  comment=FALSE
  
  fid=file(filename,'r')
  vmodini=rep(NaN,6)
  vmod=NULL
  
  parmcnt=1
  namecnt=1
  linecnt=1
  done=0
  while( done==0){
    line=readLines(fid,1)
    if( length(line)==0 ){
      done=1
    }else if(comment==TRUE & !('/' %in% strsplit(line,'*',fixed=TRUE)[[1]][2:length(strsplit(line,'*',fixed=TRUE)[[1]])]) ){
    }else if(comment==TRUE){
      comment=FALSE 
    }else{
      ## strip /**/ comments
      if(length(strsplit(line,'/*',fixed=TRUE)[[1]])>1){
        comment=TRUE # ignore lines until */ is reached
      }
      line=strsplit(line,'/*',fixed=TRUE)[[1]][1]
      ## strip // and # comments
      line=strsplit(line,'//',fixed=TRUE)[[1]][1]
      line=strsplit(line,'#',fixed=TRUE)[[1]][1]
      
      if(nchar(line)>0 & !is.na(line)){
        vec = strsplit(line,' ')[[1]]
        if(is.na(sum(as.numeric(vec[nzchar(vec)])))){
                                        # line includes letters				
          nws=paste(strsplit(line,' ')[[1]],collapse='') 
          if(substr(nws,1,1)=='!'){
            
            vec = strsplit(line,' ')[[1]]
            keyword = tolower(vec[which(nzchar(vec))[1]])
            par = vec[which(nzchar(vec))[2]]
            if(keyword=='!radius'){
              pradius=as.numeric(par)
            }else if(keyword=='!year'){
              pubyear=as.numeric(par)
            }else if(keyword=='!name'){
              mname=par
            }else{
              if(verbose){
                print(paste('ReadND: unrecognized keyword',keyword,'ignored.'))
              }
            } 
          }else{ 
            vec = strsplit(line,' ')[[1]]
            w = which(nzchar(vec))
            name = paste(vec[min(w):max(w)],collapse = ' ')
            dz[namecnt]=vmod[parmcnt-1,1]
            dname=c(dname,name)
            if(verbose){
              print(paste('ReadND:',name,'discontinuity at',dz[namecnt],'km'))
            }
            namecnt=namecnt+1
          }
        }else{ 
          
          pcnt=1
          done2=0
          vmod=rbind(vmod,vmodini,deparse.level=0)
          
          vec = strsplit(line,' ')[[1]]
          vmod[parmcnt,1:sum(nzchar(vec))] = as.numeric(vec[nzchar(vec)])
          parmcnt=parmcnt+1
        } 
      } 
      linecnt=linecnt+1
    }
  } 
  
  close(fid)
  
##### replace all -1 by NaN
  indies=which(vmod==-1)
  vmod[indies]=NaN
  indies=which(dz==-1)
  dz[indies]=-1
  
  
  if( is.null(pradius) ){
    pradius=max(vmod[,1])
  }
  
  
  
  model=list()
  model$z=vmod[,1]
  model$vp=vmod[,2]
  model$vs=vmod[,3]
  model$rho=vmod[,4]
  model$qp=vmod[,5]
  model$qs=vmod[,6]
  model$name=mname
  model$rp=pradius
  model$year=pubyear
  
  model$conr=NaN
  model$moho=NaN
  model$d410=NaN
  model$d520=NaN
  model$d660=NaN
  model$cmb=NaN
  model$icb=NaN
  model$dz=NULL
  model$dname=NULL
  anz=length(dz) 
  cnt=1 
  
  for( indy in 1:anz){
    vec = strsplit(dname[indy],' ')[[1]]
    w = which(nzchar(vec))
    p=tolower(paste(vec[min(w):max(w)],collapse=' '))
    if(p == 'conrad'){
      model$conr=dz[indy]
    }else if(p %in% c('moho','mantle')){
      model$moho=dz[indy]
    }else if(p %in% c('olivine alpha beta','transition zone')){
      model$d410=dz[indy]
    }else if(p == 'olivine beta gamma'){
      model$d520=dz[indy]
    }else if(p %in% c('olivine gamma perovskite','lower mantle')){
      model$d660=dz[indy]
    }else if(p %in% c('cmb','outer core','outer-core')){
      model$cmb=dz[indy]
    }else if(p %in% c('icb','inner core','inner-core')){
      model$icb=dz[indy]
    }else{
      if(verbose){
        print(paste('ReadND: non-standard discontinuity',p))
      }
      model$dz=c(model$dz, dz[indy])
      model$dname=c(model$dname,p)
    }
  } 
  
  if(verbose){
    print(paste('ReadND: planet radius is', pradius))
    print(paste('ReadND:', namecnt-1, 'discontinuity names.'))
    print(paste('ReadND:',parmcnt-1, 'parameter sets.'))
    print(paste('ReadND:',linecnt-1, 'lines total.'))
  }
  return(model)
}

