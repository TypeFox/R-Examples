#' @keywords internal
akestN<-function(N,s,
          K=max(s),
          verbose=FALSE, 
          n=NULL,
	  maxit=5000, tol=0.00000001,
          return.all=FALSE){
  if(length(s)>N){print("Error: The population counts should not be exceeded by the sample counts.")}
# sizes<-sort(unique(s))
  sizes<-1:K
  s[s>K] <- K
  if(is.null(n)){n=tabulate(s,nbins=K)}
  enroot<-function(llam,N,nk,sizes){
    sum(nk/(1-exp(-exp(llam)*sizes)))-N
  }
  out <- uniroot(f=enroot, interval=c(-20,5),
     N=N, nk=n, sizes=sizes,
#    f.lower=Inf, f.upper=sum(nk)-N,
     maxiter=maxit,tol=tol)
  out$lambda<-exp(out$root)
  out$NNhat<-n/(1-exp(-out$lambda*sizes))
  if(return.all){
    out
  }else{
    out$NNhat
  }
}
#' @keywords internal
sppsestN<-function(N,s,
          K=max(s),
          verbose=FALSE, 
          n=NULL,
	  maxit=10, M=10000, trace=FALSE,
          return.all=FALSE){
  if(length(s)>N){print("Error: The population counts should not be exceeded by the sample counts.")}
# sizes<-sort(unique(s))
  s[s>K] <- K
  if(is.null(n)){n=tabulate(s,nbins=K)}
  out <- getestCstacked(samp=s,n=N,nit=maxit,nsampsamp=M)
  NNhat <- n
  NNhat[n!=0] <- out$NNhat
# out$NNhat.raw <- out$NNhat
  out$NNhat <- NNhat
  if(return.all){
    out
  }else{
    out$NNhat
  }
}
#' @keywords internal
roundstoc<-function(vec){# takes a vector and makes it integers, keeping the total the same
	#vec<-c(1.2,1.3,2.5)
	target<-sum(vec)
	temp<-floor(vec)
	needed<-target-sum(temp)
	away<-vec-temp
	while(needed>.5){
		toget<-sample(c(1:length(vec)),size=1,prob=away)
		temp[toget]<-temp[toget]+1
		away<-vec-temp
		away[away<0]<-0
		needed<-needed-1
		}
	temp
	}

#' @keywords internal
probtodist<-function(classes,nums,prob,n){
  #classes: size value of each class 
  #nums: number of members of that class 
  #prob: Probabilities of selection for a member of that class
	props<-rep(0,length(classes))
	for(k in 1:length(classes)){
	  props[k]<-nums[k]/prob[k]
	}
	props<-props/sum(props)
	if(any(is.na(props))){print("Error in proptodist")}	
  #classes: size value of each class: 
  #props: proportion of the population in that class
	list(classes=classes, props=props)	
	}


#' @keywords internal
getincl<-function(classes,props,n,nsamp,nsampsamp){
	#print(paste("n starting getincl is",n))
# 2a) Estimate set of degrees for all nodes
props2<-props/sum(props)
nbyclass<-round(n*props2)
#artificially inflate each value to at least 1 node
nbyclass[nbyclass==0]<-1
offby<-n-sum(nbyclass)
#print(paste("offby=",offby))
if(is.na(offby)){print("Error in getincl: offby"); 
	#debug()
	}
#print(offby)	##
if(offby>0){
	for(ii in 1:offby){
		dif<-props2-nbyclass/n
		dif[dif<0]<-0
		tochange<-sample(1:length(dif),1,prob=dif)
		nbyclass[tochange]<-nbyclass[tochange]+1
#		nbyclass[which.max(temp2)]<-nbyclass[which.max(temp2)]+1
		}
}else{if(offby<0){
	for(ii in 1:abs(offby)){		
#	  	temp2<-props-nbyclass/n
		dif<-nbyclass/n - props2
		#print(dif) ##
		dif[dif<0]<-0
		#print(dif) ##
		dif[nbyclass==1]<-0 #don't get rid of the last of any class
        if(sum(dif==0)){
        	dif<-1-(props2-nbyclass/n)
        	#print("In dif=0 loop")
        	dif[nbyclass==1]<-0
        	#print(cbind(dif,nbyclass))
        	}	
		#print(dif) ##
		#print(cbind(dif,nbyclass))
		tochange<-sample(1:length(dif),1,prob=dif)
		nbyclass[tochange]<-nbyclass[tochange]-1
#		nbyclass[which.min(temp2)]<-nbyclass[which.min(temp2)]		
        }
	}}
	#print(min(nbyclass))
popclass<-rep(classes,times=nbyclass)	

# based on estimate of degrees, estimate true inclusion probs by class
  timesdeg<-rep(0,n)
  timesincldeg<-rep(0,n)
    includes<-rep(0,n)
    for(i in 1:nsampsamp){
  	  aaa<-sample(c(1:n),size=(nsamp),replace=FALSE, 
  	    prob=popclass)
  	    #print(paste("n in second loop is",n))
  	    #print(paste("length of popclass in second loop is",length(popclass)))
  	  includes[aaa]<-includes[aaa]+1
  	  }
      # keep running tally of number of nodes of each degree and   
      # number of inclusions of nodes of each degree
     for(i in 1:n){
      bbb<-popclass[i]
   	  if(bbb>0){
   	    timesdeg[bbb]<-timesdeg[bbb]+1
   	    timesincldeg[bbb]<-timesincldeg[bbb]+includes[i]
        } 
   	  }
    degvec<-classes #which(timesdeg>0)
    #artificially inflate sample so each is seen at least once
    timesincldeg[classes][timesincldeg[classes]==0]<-timesincldeg[classes][timesincldeg[classes]==0]+1
    #print(classes)
    #print(timesdeg[classes])
    pvec<-timesincldeg[classes]/timesdeg[classes]
    if(any(is.na(pvec))){print("Error: 0 probability")}
    #print(pvec)
    pvec<-pvec/sum(pvec)
    #print(pvec)
    if(any(is.na(pvec))){print("Error in incldeg: pvec")}
    list(degvec=degvec,pvec=pvec,nbyclass=nbyclass)
}

#' @keywords internal
getest<-function(samp,n,nit=5, nsampsamp=1000, trace=FALSE,remember=FALSE){
  classes<-sort(unique(samp))
  nums<-table(samp)
  nsamp<-length(samp)
  prob<-classes/sum(samp)
  memory<-matrix(0,nit+1,length(classes))

  temp<-probtodist(classes,nums,prob,n)
  newprobs<-getincl(temp$classes,temp$props,n, nsamp, nsampsamp)

  for(i in 1:nit){
    temp<-probtodist(newprobs$degvec,nums,newprobs$pvec,n)
    newprobs<-getincl(temp$classes,temp$props, n, nsamp, nsampsamp=nsampsamp)
  }
  list(probs=newprobs$pvec,NNhat=n*temp$props,sizes=temp$classes)
}

#' @keywords internal
getincl.raw<-function(size,
                  Nk=sort(unique(size)),
		  n=500,
                  M=1000,
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    K <- length(Nk)
    pop <- match(size,Nk)-1
    if(!is.null(seed))  set.seed(as.integer(seed))
    size <- size/sum(size)
    Cret <- .C("getinclC",
              N=as.integer(length(size)),
              pop=as.integer(pop),
              size=as.double(size),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(M),
              Nk=as.integer(Nk), PACKAGE="sspse")
#             verbose=as.integer(verbose))
    Cret
}
