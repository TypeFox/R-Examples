OptimizeDist <- function(alphalimit,deltalimit,phase,h,imodel){
  
### Define internal function:
  ConvAng2pWrap=function(phase,h,angle,vp,vs,rp){
    
    if( length(vp)==2){ 
      if( angle<=90){
        whichvelocity=2
      }else{
        whichvelocity=1
      } 
    }else{
      whichvelocity=1
    } 
    
    p=ConvAng2p(phase,h,angle,model=NULL,vp=vp[whichvelocity],vs=vs[whichvelocity],rp=rp)
    return(p)
  }
  
  QualifyExtremum=function(a,d){ # is the center of 3 points a local min, max, or neither
    sorter=order(a)
    a=a[sorter]
    d=d[sorter]
    
    if( (d[2] < d[1]) && (d[2] < d[3])){
      qualifier=-1 # local minimum
    }else{
      if( (d[2] > d[1]) && (d[2] > d[3])){
        qualifier=1 # local maximum
      }else{
        qualifier=0 # neither
      }
    }
    return(qualifier)
  }
  
### Done defining internal functions
  
  tt=proc.time()
  eps = .Machine$double.eps
  extremalpha=NULL
  extremp=NULL
  extremdelta=NULL

  sorter=order(alphalimit)
  alphalimit=alphalimit[sorter]
  deltalimit=deltalimit[sorter]

  focus=InterpModel(imodel,h,'simple')
  indy=which(focus$z==h)
  vp=focus$vp[indy]
  vs=focus$vs[indy]
  identepsilon=1e-6 # changed from 1e-8 because sometimes diff(deltalimit)==0 while diff(alphalimit) > 1e-8
  
  if( (diff(alphalimit) < identepsilon) && (diff(deltalimit) < identepsilon) ){
    extremalpha=alphalimit[1]
    extremdelta=deltalimit[1]
    extremp=ConvAng2pWrap(phase, h, extremalpha, vp, vs, imodel$rp)
    return(list(extremalpha=extremalpha,extremp=extremp,extremdelta=extremdelta))
  } 
  phi=(1+sqrt(5))/2 
  goldrat=phi-1 
  gradientstep=0.00001 
  goldepsilon=eps^.25
  
  
  goldmaxcnt=ceiling((log(alphalimit[length(alphalimit)]-log(alphalimit[1]))-log(goldepsilon))/log(phi))
  goldmaxcnt=goldmaxcnt+1 
  alphalist=alphalimit 
  deltalist=deltalimit 
  alphalen=length(alphalist)
  rayplist=alphalist * 0
  for(i in 1:length(alphalist)){ # This for loop was added 5/16/11--previously it just had rayplist = CongAng2pWrap(phase, h, alphalist, vp, vs, imodel$rp), which caused a warning.
    rayplist[i] = ConvAng2pWrap(phase, h, alphalist[i], vp, vs, imodel$rp)
  }
  a=NULL
  a[1]=alphalimit[1] 
  a[3]=alphalimit[2] 
  a[2]=a[1]+goldrat*(a[3]-a[1]) 
  p=NULL
  p[1]=rayplist[1]
  p[3]=rayplist[2]
  p[2]=ConvAng2pWrap(phase,h,a[2],vp,vs,imodel$rp)
  d=NULL
  d[1]=deltalimit[1]
  d[3]=deltalimit[2]
  d[2]=FindDist4p(phase,h,imodel,p=p[2],takeoff=a[2])[[1]]
  extremumfound=0 
  done=FALSE
  cnt=0
  while(!done){
    
    
    alphalist=c(alphalist, a[2])
    deltalist=c(deltalist, d[2])
    rayplist=c(rayplist, p[2])
    
    if(sum(is.na(c(p,a,d)))){return(list(NaN,NaN,NaN))}
    
    qualifier=QualifyExtremum(a,d)

    if(qualifier == -1){
      searchwhat=qualifier
      done=TRUE
    }else if(qualifier == 1){
      searchwhat=qualifier
      done=TRUE
    }else if(qualifier == 0){
      searchwhat=0
      if (cnt==0){
        a[2]=a[3]-goldrat*abs(a[3]-a[1])
        p[2]=ConvAng2pWrap(phase,h,a[2],vp,vs,imodel$rp)
        d[2]=FindDist4p(phase,h,imodel,p=p[2],takeoff=a[2])[[1]]
        
        alphalist=c(alphalist,a[2])
        deltalist=c(deltalist,d[2])
        rayplist=c(rayplist,p[2])
      }else{
                
        intervallength=sqrt((a[3]-a[1])^2+(d[3]-d[1])^2) # JFA can't imagine why one would use the 2-D distance here, but it shouldn't matter and if it aint broke, don't fix it...5/18/11
        rightofa=a[1]+(gradientstep*intervallength) 
        leftofc=a[3]-(gradientstep*intervallength) 
        
        prightofa=ConvAng2pWrap(phase,h,rightofa,vp,vs,imodel$rp) 
        pleftofc=ConvAng2pWrap(phase,h,leftofc,vp,vs,imodel$rp) 
        aa=FindDist4p(phase,h,imodel,p=c(prightofa,pleftofc),takeoff=c(rightofa,leftofc))$dist
        frightofa=aa[1] 
        fleftofc=aa[2] 
        if(is.na(frightofa)|is.na(fleftofc)){return(list(NaN,NaN,NaN))}
        grada=(frightofa-d[1])/(gradientstep*intervallength) 
        gradc=(d[3]-fleftofc)/(gradientstep*intervallength) 
        
        alphalist=c(alphalist, rightofa, leftofc)
        deltalist=c(deltalist, frightofa,fleftofc)
        rayplist=c(rayplist, prightofa, pleftofc)
        fdiscr=sign(d[1]-d[3]) 
        
        if(fdiscr == -1){
          if( grada>0){
            if( gradc>=0 & !is.na(gradc)){
              searchwhat = 1
              extremalpha=a[3]
              extremp=p[3]
              extremdist=d[3]
              done=1
              extremumfound=1
            }else{
              searchwhat=+1
              a[2]=rightofa
              p[2]=prightofa
              d[2]=frightofa
            } 
          }else{
            
            if( gradc>0 & !is.na(gradc)){
              searchwhat=-1
              a[2]=rightofa
              p[2]=prightofa
              d[2]=frightofa
              done=1
            }else{
              return(list(extremalpha=NaN,extremdelta=NaN,extremp=NaN))
            } 
          } 
          
          
        }else if(fdiscr == 0){
          
          if(sign(grada)== -1){
            searchwhat=-1
            a[2]=rightofa
            p[2]=prightofa
            d[2]=frightofa
            done=1
          }else if(sign(grada)==0){
            stop('OptimizeDist: insufficient information!')
          }else if(sign(grada)==1){
            searchwhat=+1
            a[2]=rightofa
            p[2]=prightofa
            d[2]=frightofa
            done=1
          }else{
            stop('OptimizeDist: unexpected gradient sign (1)!')
          } 
        }else if(fdiscr==1){
          
          if(sign(grada)==-1){
            if( gradc>=0 && !is.na(gradc)){
              searchwhat=-1
              a[2]=leftofc
              p[2]=pleftofc
              d[2]=fleftofc
              done=1
            }else{
              searchwhat=-1
              extremalpha=a[3]
              extremp=p[3]
              extremdist=d[3]
              done=1
              extremumfound=1
            } 
          }else if(!sign(grada)){
            
            searchwhat=-1
            extremalpha=a[1]
            extremp=p[1]
            extremdist=d[1]
            done=1
            extremumfound=1
          }else if(sign(grada)== 1){
            searchwhat= 1
            a[2]=rightofa
            p[2]=prightofa
            d[2]=frightofa
            done=1
          }else{
            stop('OptimizeDist: unexpected gradient sign (2)!')
          } 
        }else{
          warning(paste('OptimizeDist: unexpected exception. Gradients: A: ', (grada), ', C: ',(gradc),' (Delta=', (deltalimit[1]), '...', (deltalimit[2]), ')'))
        } 
      }
    } 
    
    cnt=cnt+1
    
    done=(cnt>2) | done
    extremumfound=(cnt>2) | extremumfound
  } 
  sorter=order(a)
  a=a[sorter]
  p=p[sorter]
  d=d[sorter]
  if( extremumfound==0){
    done=FALSE
    cnt=0
    while(!done){
      intlen=abs(diff(a))
      if( sum(intlen)>goldepsilon){
        if( intlen[1]>=intlen[2]){
          newa=a[1]+goldrat*intlen[1]
        }else{
          newa=a[3]-goldrat*intlen[2]
        } 
        
        newp=ConvAng2pWrap(phase,h,newa,vp,vs,imodel$rp) 
        newd=FindDist4p(phase,h,imodel,p=newp,takeoff=newa)[[1]] 
        
        alphalist=c(alphalist, newa)
        rayplist=c(rayplist, newp)
        deltalist=c(deltalist,newd)
        
        a=c(a, newa)
        p=c(p, newp)
        d=c(d, newd)

        sorter=order(a)
        a=a[sorter]
        p=p[sorter]
        d=d[sorter]
        
        if(searchwhat== -1){
          extremindy=which(min(d)==d)[1]
        }else if(searchwhat == 1){
          extremindy=which(max(d)==d)[1]
        }else{
          stop(paste('OptimizeDist: unexpected SEARCHWHAT: ', (searchwhat)))
        } 
        if(extremindy==1){
          newpoints=1:3
        }else if(extremindy==4){
          newpoints=2:4
        }else{
          newpoints=extremindy + (-1):1
        } 
        a=a[newpoints]
        p=p[newpoints]
        d=d[newpoints]
        
      }else{
        done=TRUE
      } 
      
      cnt=cnt+1
      if( cnt>goldmaxcnt){
        done=TRUE
        print('OptimizeDist: iteration exhaust!')
      } 
    } 
  } 
  if(searchwhat==-1){
    extremdelta=min(deltalist)
    indy=(deltalist==extremdelta)
  }else if(searchwhat==1){
    extremdelta=max(deltalist)
    indy=(deltalist==extremdelta)
  }else if(searchwhat==0){
    indy=NULL
  }else{
    stop('OptimizeDist: unexpected SEARCHWHAT value.')
  } 
  return(list(extremalpha=(alphalist[indy])[1],extremp=(rayplist[indy])[1],extremdelta=extremdelta))
}

