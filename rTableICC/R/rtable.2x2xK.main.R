rtable.2x2xK.main <-
function(p,sampling,N,K,lambda,print.raw){
  if (sampling=="Product"){
    if ((length(N)==1) | (is.finite(N)==FALSE)){ 
      stop("Total number of observations under each center should be entered as a finite Kx1 vector under product multinomial samlping plan.")      
    }
    N=abs(round(N))  
    center.margins=N
    K=length(center.margins)
  }else if (sampling=="Multinomial"){ 
    K=dim(p)[3]
  }else if ((sampling=="Poisson") & (is.null(K)==TRUE)){
    K=dim(lambda)[3]
  }
  
  rTable=array(0,dim=c(2,2,K))  

  if (sampling=="Product"){
    N=sum(center.margins)    
    rTable.raw=array(0,dim=c(N,3))  
    gen=array(0,4)        
    for (i in 1:K){
      p.row=as.vector(t(p[,,K]/(center.margins[K]/N)))
      gen=rmultinom(1, center.margins[i], p.row)  
      say=0   
      bsl=1 
      for (l in 1:2){
        for (r in 1:2){
          say=say+1
          rTable[l,r,i]=gen[say]
          if (rTable[l,r,i]>0){
            bts=(bsl+rTable[l,r,i]-1)
            rTable.raw[bsl:bts,1]=l
            rTable.raw[bsl:bts,2]=r
            rTable.raw[bsl:bts,3]=i
            bsl=bsl+rTable[l,r,i]          
          }
        }
      }
    }  
  } else if (sampling=="Multinomial"){
    
      N=abs(round(N)) 
      pp=array(0,2*2*K)
      say=0
      for (i in 1:K){
        for (l in 1:2){
          for (r in 1:2){
            say=say+1
            pp[say]=p[l,r,i]
          }
        }
      }
      rTable.raw=array(0,dim=c(N,3)) 
      gen=rmultinom(1,N,pp)  
      say=0
      bsl=1 
      for (i in 1:K){
        for (l in 1:2){
          for (r in 1:2){
            say=say+1
            rTable[l,r,i]=gen[say]
            if (rTable[l,r,i]>0){
              bts=(bsl+rTable[l,r,i]-1)
              rTable.raw[bsl:bts,]=l
              rTable.raw[bsl:bts,2]=r
              rTable.raw[bsl:bts,3]=i
              bsl=bsl+rTable[l,r,i]          
            }
          }
        }
      }    
  } else if (sampling=="Poisson"){ 
    if (length(lambda)>1){
      pp=array(0,2*2*K)
      say=0
      for (i in 1:K){
        for (l in 1:2){
          for (r in 1:2){
            say=say+1
            pp[say]=lambda[l,r,i]
          }
        }
      }
      gen=rpois(2*2*K,pp)  
    } else {
      gen=rpois(2*2*K,lambda)
    } 
    N=sum(gen)          
    rTable.raw=array(0,dim=c(N,3))      
    say=0
    bsl=1 
    for (i in 1:K){
      for (l in 1:2){
        for (r in 1:2){
          say=say+1
          rTable[l,r,i]=gen[say]
          if (rTable[l,r,i]>0){
            bts=(bsl+rTable[l,r,i]-1)
            rTable.raw[bsl:bts,1]=l
            rTable.raw[bsl:bts,2]=r
            rTable.raw[bsl:bts,3]=i
            bsl=bsl+rTable[l,r,i]          
          }
        }
      }
    }    
  } 

  list(rTable=rTable,rTable.raw=rTable.raw,N=N,sampling=sampling,K=K,ICC=FALSE,structure="2x2xK",print.raw=print.raw)
    
}
