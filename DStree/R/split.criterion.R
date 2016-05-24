temp1S <- function(y, wt, parms) {
  n <- length(parms[[1]])
  dat.haz <- hazard(y[,1],y[,2],parms[[1]],wt)
  lab <- rep(NA,2*n+1)
  dev <- -2*(likelihood(dat.haz))
  
  S.1 <- rev(dat.haz$S[dat.haz$S>=0.5])[1]
  S.2 <- dat.haz$S[dat.haz$S<0.5][1]
  
  if(is.na(S.1)){
  S.1 <- 1
  m <- 0
  }else{
  ms <- which(dat.haz$S==S.1)
  m <- ms[length(ms)]
  }
  
  
  if(is.na(S.2)) {
    S<-NA
  } else {
  S <-round(m+(S.1-0.5)/S.1-S.2,1)
  }
  
  
  lab[1] <- S
  lab[2:(n+1)] <- dat.haz$h
  lab[(n+2):(2*n+1)] <-dat.haz$S
  list(label=lab, deviance=dev)
}

temp2S <- function(y, wt, x, parms, continuous) {
   if (continuous) {
    n <- length(y[,1])
    duplicate <- !duplicated(x)
    uniq.ind <- seq_along(x)[duplicate]
    k<-length(uniq.ind)
    
    temp.llik<-temp.rlik<-rep(NA,k-1)
    for(i in 2:k){
      surv.l <- lapply(uniq.ind[i]-1,hazardl,y2=y[,2],lev=parms[[1]],y1=y[,1],wt=wt)[[1]]
      temp.llik[i-1]<-lik(surv.l$S,surv.l$pi,surv.l$n.cens,surv.l$n.uncens)
      surv.r <- lapply(uniq.ind[i],hazardr,y2=y[,2],lev=parms[[1]],y1=y[,1],wt=wt)[[1]]
      temp.rlik[i-1]<-lik(surv.r$S,surv.r$pi,surv.r$n.cens,surv.r$n.uncens)
    }
    llik<-n+temp.llik/n
    rlik<-n+temp.rlik/n
    
    goodness <-rep(0,n-1)
    goodness[uniq.ind[-1]-1] <- llik+rlik
    list(goodness=goodness , direction=rep(1,n-1))
    
  }
  else {
    ux <- sort(unique(x))
    m <- length(x)
    y.nlev<-length(unique(y[,2]))
    M=cbind(y[,1],x,y[,2],wt)
    
   
    means <- rep(0,length(ux))
    for(i in 1:length(ux)){
      
      dat.sub <- matrix(M[M[,2]==i,c(1,3,4)],ncol=3)
      dat.haz <- hazard(dat.sub[,1],dat.sub[,2],parms[[1]],dat.sub[,3])
      means[i] <- likelihood(dat.haz)
      
    }
    
    ord <- order(means)
    n <- length(ord)
    llik <- rlik <- rep(0,n-1)
    
    for(i in 1:(n-1)){
      dat.2 <- matrix(M[M[,2] %in% ux[ord[1:i]],c(1,3,4)],ncol=3)
      dat.3 <- matrix(M[M[,2] %in% ux[ord[(i+1):n]],c(1,3,4)],ncol=3)
      
      dat.haz.2 <- hazard(dat.2[,1],dat.2[,2],parms[[1]],dat.2[,3])
      dat.haz.3 <- hazard(dat.3[,1],dat.3[,2],parms[[1]],dat.3[,3])
      llik[i] <- m+likelihood(dat.haz.2)/m
      rlik[i] <- m+likelihood(dat.haz.3)/m
      
    }
    
    list(goodness= llik + rlik,
         direction = ux[ord])
    
  }
}



temp3S <- function(y, offset, parms, wt) {
  
  unique.y <-sort(unique(y[,1]))
  numresp=2*length(unique.y)+1
  sfun <- function(yval, dev, wt, ylevel, digits ) {
    paste("Median Surv=", format(signif(yval, digits)),
          ", Deviance=" , format(signif(dev, digits)),
          ", wt=" , format(signif(wt, digits)),
          # ", ylev=" , format(str(ylevel)),
          sep='')}
  environment(sfun) <- .GlobalEnv
  
  list(y=y, parms=list(unique.y), numresp=numresp, numy=2,summary=sfun)
}

alist <- list(eval=temp1S, split=temp2S, init=temp3S)
