rtable.RxC.main <-
function(p,row.margins,col.margins,sampling,N,lambda,print.raw){
  
  if (sampling=="Product"){
 
    R=nrow(p)
    C=ncol(p)
    rTable=array(0,dim=c(R,C))  
    
    if (length(row.margins)>1){
      N=sum(row.margins)
      prob.margins=row.margins/N
      if (round(apply(p,1,sum),10)!=round(prob.margins,10)){
        stop("Mismatch between cell probabilities and margin probabilities!")
      }     
      p=sweep(p, 1, prob.margins, "/")  
      rTable.raw=array(0,dim=c(N,2))  
      bsl=1
      for (i in 1:R){
        rTable[i,1:C]=rmultinom(1, row.margins[i], t(p)[(((i-1)*C)+1):(i*C)]) 
        for (j in 1:C){
          if (rTable[i,j]>0){
            rTable.raw[bsl:(bsl+rTable[i,j]-1),1]=i
            rTable.raw[bsl:(bsl+rTable[i,j]-1),2]=j
            bsl=bsl+rTable[i,j]
          }
        }        
      }
    }else if (length(col.margins)>1){
      N=sum(col.margins)
      prob.margins=col.margins/N
      if (round(apply(p,1,sum),10)!=round(prob.margins,10)){
        stop("Mismatch between cell probabilities and margin probabilities!")
      }      
      p=sweep(p, 2, prob.margins, "/")
      rTable.raw=array(0,dim=c(N,2))  
      bsl=1
      for (i in 1:C){
        rTable[1:R,i]=rmultinom(1, col.margins[i], p[(((i-1)*R)+1):(i*R)])  
        for (j in 1:R){
          if (rTable[j,i]>0){
            rTable.raw[bsl:(bsl+rTable[j,i]-1),1]=j
            rTable.raw[bsl:(bsl+rTable[j,i]-1),2]=i
            bsl=bsl+rTable[j,i]
          }
        }        
      }
    }    
  } else if (sampling=="Multinomial"){     
    if ((length(N)!=1) | (is.finite(N)==FALSE)){ 
      stop("Total number of observation should be entered as a scalar under multinomial samlping plan.")      
    }
    N=abs(round(N))
    R=nrow(p)
    C=ncol(p)
    rTable=array(0,dim=c(R,C)) 
    rTable.raw=array(0,dim=c(N,2))
    p=rmultinom(1, N, p) 
    bsl=1
    say=0
    for (i in 1:R){         
      for (j in 1:C){
        say=say+1
        rTable[i,j]=p[say]          
          if (rTable[i,j]>0){
            rTable.raw[bsl:(bsl+rTable[i,j]-1),1]=i
            rTable.raw[bsl:(bsl+rTable[i,j]-1),2]=j
            bsl=bsl+rTable[i,j]
          }
      }        
    }      
  }else if (sampling=="Poisson"){  
    R=nrow(lambda)
    C=ncol(lambda)
    rTable=array(0,dim=c(R,C))  
    p=rpois(R*C,t(lambda))  
    N=sum(p)
    rTable.raw=array(0,dim=c(N,2))
    bsl=1
    say=0
    for (i in 1:R){         
      for (j in 1:C){
        say=say+1
        rTable[i,j]=p[say] 
        if (rTable[i,j]>0){
          rTable.raw[bsl:(bsl+rTable[i,j]-1),1]=i
          rTable.raw[bsl:(bsl+rTable[i,j]-1),2]=j
          bsl=bsl+rTable[i,j]
        }
      }        
    }  
  }
  list(rTable=rTable,rTable.raw=rTable.raw,N=N,sampling=sampling,R=R,C=C,ICC=FALSE,structure="RxC",print.raw=print.raw)
}
