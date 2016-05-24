

#############################################
#####  Prepare data for log-likelihood  #####
#############################################

coxprep=function(x, y){
  
  N0=nrow(x)
  oi=order(y[, "status"], decreasing=TRUE)
  x=x[oi, ];y=y[oi, ]
  oi=order(y[, "time"])
  x=x[oi, ];y=y[oi, ]
  
  ## remove the first censored cases
  i1=which(y[, "status"]==1);mi1=min(i1)-1
  if (mi1!=0) {
    x=x[-c(1:mi1), ];y=y[-c(1:mi1), ]
  }
  ty=y[, "time"];tevent=y[, "status"]
  N=nrow(x);n1=sum(y[, "status"])
  
  dty=duplicated(ty) # ties
  
  ### for calculation of log-likelihood
  if (any(dty)) {
    tevent0=tevent
    tevent0[which(dty)]=0
    
    ievent=cumsum(tevent0);loc1=which(tevent0==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=tapply(tevent==1, ievent, sum)
  } else {
    ievent=cumsum(tevent);loc1=which(tevent==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=rep(1, n)
  }
  
  return(list(x=x, N0=N0, tevent=tevent, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n))
}



###############
###  local  ###
###############

locoxprep=function(x, y, w, w0, h){
  wi=(w-w0)/h;indexi=(abs(wi)<=1)
  
  if (sum(indexi)>=2) {
    if (any(y[indexi, "status"]==1)) {
      N0=sum(indexi) # number of effective data
      Kh=ifelse(indexi, (1-wi^2)*3/4, 0)/h # Epanechnikov
      x=x[which(indexi), ];y=y[which(indexi), ];Kh=Kh[which(indexi)]
      
      oi=order(y[, "status"], decreasing=TRUE)
      x=x[oi, ];y=y[oi, ];Kh=Kh[oi]
      oi=order(y[, "time"])
      x=x[oi, ];y=y[oi, ];Kh=Kh[oi]
      
      ## remove the first censored cases
      i1=which(y[, "status"]==1);mi1=min(i1)-1
      if (mi1!=0) {
        x=x[-c(1:mi1), ];y=y[-c(1:mi1), ];Kh=Kh[-c(1:mi1)]
      }
      
      if (length(Kh)>1) {
        ty=y[, "time"];tevent=y[, "status"];N=nrow(x)
        
        dty=duplicated(ty) # ties
        ### for calculation of log-likelihood
        if (any(dty)) {
          tevent0=tevent
          tevent0[which(dty)]=0
          
          ievent=cumsum(tevent0);loc1=which(tevent0==1)
          nevent=table(ievent);n=length(unique(ievent))
          nevent1=tapply(tevent==1, ievent, sum)
          Kh1=tapply(tevent*Kh, ievent, sum)
        } else {
          ievent=cumsum(tevent);loc1=which(tevent==1)
          nevent=table(ievent);n=length(unique(ievent))
          nevent1=rep(1, n);Kh1=Kh[loc1]
        }
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
  
  return(list(x=x, N0=N0, tevent=tevent, Kh=Kh, Kh1=Kh1, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n))
}




