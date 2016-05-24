#' @keywords internal
getestC<-function(samp,n,nit=5, nsampsamp=1000){
# samp: sample
# n: population size
  classes<-sort(unique(samp))
  nums<-table(samp)
  nsamp<-length(samp)
  prob<-classes/sum(samp)

  #temp$classes: size value of each class 
  #temp$nums: number of members of that class 
  #temp$prob: Probabilities of selection for a member of that class
  temp<-probtodist(classes=classes,nums=nums,prob=prob,n=n)
  #temp$classes: size value of each class 
  #temp$props: proportion of the population in that class
  newprobs<-getinclC(classes=temp$classes,props=temp$props,n=n, nsamp=nsamp,
                     nsampsamp=nsampsamp)
  #temp$degvec: size value of each class
  #temp$pvec: Probabilities of selection for a member of that class

  for(i in 1:nit){
    temp<-probtodist(classes=newprobs$degvec,nums=nums,prob=newprobs$pvec,n=n)
    newprobs<-getinclC(classes=temp$classes,props=temp$props,n=n, nsamp=nsamp,
                       nsampsamp=nsampsamp)
  }
  #classes: size value of each class 
  #props: proportion of the population in that class
  #probs: Probabilities of selection for a member of that class
  list(probs=newprobs$pvec,props=temp$props,classes=temp$classes)
}
#' @keywords internal
getinclC<-function(classes,props,n,nsamp,nsampsamp){
  props2<-props/sum(props)
  nbyclass<-round(n*props2)
  #artificially inflate each value to at least 1 node
  nbyclass[nbyclass==0]<-1
  offby<-n-sum(nbyclass)
  if(is.na(offby)){print("Error in getincl: offby")}
  if(offby>0){
   for(ii in 1:offby){
    dif<-props2-nbyclass/n
    dif[dif<0]<-0
    tochange<-sample(1:length(dif),1,prob=dif)
    nbyclass[tochange]<-nbyclass[tochange]+1
   }
  }else{if(offby<0){
   for(ii in 1:abs(offby)){		
    dif<-nbyclass/n - props2
    dif[dif<0]<-0
    dif[nbyclass==1]<-0 #don't get rid of the last of any class
    if(sum(dif==0)){
     dif<-1-(props2-nbyclass/n)
     dif[nbyclass==1]<-0
    }	
    tochange<-sample(1:length(dif),1,prob=dif)
    nbyclass[tochange]<-nbyclass[tochange]-1
   }
  }}
  popclass<-rep(classes,times=nbyclass)	
# based on estimate of degrees, estimate true inclusion probs by class
  pop <- match(popclass,classes)-1
  popclass <- popclass/sum(popclass)
  Cret <- .C("getinclC",
              N=as.integer(length(popclass)),
              pop=as.integer(pop),
              size=as.double(popclass),
              K=as.integer(length(classes)),
              n=as.integer(nsamp),
              samplesize=as.integer(nsampsamp),
              Nk=as.integer(classes), PACKAGE="sspse")
#             verbose=as.integer(0))
# list(degvec=classes,pvec=Cret$Nk/(nbyclass*nsampsamp),nbyclass=nbyclass)
  Ninf <- Cret$Nk
  Ninf[Ninf==0] <- 1
  Ninf <- Ninf/nbyclass
  Ninf <- Ninf/sum(Ninf)
  list(degvec=classes,pvec=Ninf,nbyclass=nbyclass)
}

#' @keywords internal
getestCstacked<-function(samp,n,nit=5, nsampsamp=1000){
# samp: sample
# n: population size
  classes<-sort(unique(samp))
  nums<-table(samp)
  nsamp<-length(samp)
  prob<-classes/sum(samp)

  #temp$classes: size value of each class 
  #temp$nums: number of members of that class 
  #temp$prob: Probabilities of selection for a member of that class
  temp<-probtodist(classes=classes,nums=nums,prob=prob,n=n)
  #temp$classes: size value of each class 
  #temp$props: proportion of the population in that class
  newprobs<-getinclCstacked(classes=temp$classes,props=temp$props,n=n,
		     nsamp=nsamp,
                     nsampsamp=nsampsamp)
  #temp$degvec: size value of each class
  #temp$pvec: Probabilities of selection for a member of that class

  for(i in 1:nit){
    temp<-probtodist(classes=newprobs$degvec,nums=nums,prob=newprobs$pvec,n=n)
    newprobs<-getinclCstacked(classes=temp$classes,props=temp$props,n=n,
		       nsamp=nsamp,
                       nsampsamp=nsampsamp)
  }
  #classes: size value of each class 
  #props: proportion of the population in that class
  #probs: Probabilities of selection for a member of that class
  list(probs=newprobs$pvec,props=temp$props,classes=temp$classes)
}
#' @keywords internal
getinclCstacked<-function(classes,props,n,nsamp,nsampsamp){
  props2<-props/sum(props)
  nbyclass<-round(n*props2)
  #artificially inflate each value to at least 1 node
  nbyclass[nbyclass==0]<-1
  offby<-n-sum(nbyclass)
  if(is.na(offby)){print("Error in getincl: offby")}
  if(offby>0){
   for(ii in 1:offby){
    dif<-props2-nbyclass/n
    dif[dif<0]<-0
    tochange<-sample(1:length(dif),1,prob=dif)
    nbyclass[tochange]<-nbyclass[tochange]+1
   }
  }else{if(offby<0){
   for(ii in 1:abs(offby)){		
    dif<-nbyclass/n - props2
    dif[dif<0]<-0
    dif[nbyclass==1]<-0 #don't get rid of the last of any class
    if(sum(dif==0)){
     dif<-1-(props2-nbyclass/n)
     dif[nbyclass==1]<-0
    }	
    tochange<-sample(1:length(dif),1,prob=dif)
    nbyclass[tochange]<-nbyclass[tochange]-1
   }
  }}
# based on estimate of degrees, estimate true inclusion probs by class
  probclass <- classes*nbyclass/sum(classes*nbyclass)
  Cret <- .C("getinclCstacked",
              nbyclass=as.integer(nbyclass),
              size=as.double(probclass),
              K=as.integer(length(classes)),
              n=as.integer(nsamp),
              samplesize=as.integer(nsampsamp),
              Nk=as.integer(classes), PACKAGE="sspse")
  Ninf <- Cret$Nk
  Ninf[Ninf==0] <- 1
  Ninf <- Ninf/nbyclass
  Ninf <- Ninf/sum(Ninf)
  list(degvec=classes,pvec=Ninf,nbyclass=nbyclass)
}
