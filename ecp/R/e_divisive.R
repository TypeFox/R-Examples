# Create energy environment which is used by the permutation test, and change point
# estimation procedures.


#energy = new.env(parent = emptyenv())
#energy$done_ = NULL
#energy$upper = NULL
#energy$lower = NULL

#Set matrix dimensions and initialize to the zero matrix.
setDim = function(r,c,env){
	env$done_ = matrix(0, nrow=r, ncol=c)
	#assign('done_',matrix(0,nrow=r,ncol=c),envir=energy)
	#energy$done_ = matrix(0,nrow=r,ncol=c)
}

#Set the value of individual matrix entries.
setVal = function(i,j,val,env){
	old = env$done_[i,j]
	env$done_[i,j] = val
	invisible(old)
	#tmp = energy$done_
	#tmp[i,j] = val
	#assign('done_',tmp,envir=energy)
	#energy$done_[i,j] = val
}

#Functions that are used to terminate the permutation 
#test early. As of now they have been commented out since
#their backend C++ code is causing memory problems.

#setUpper = function(x,env){
#	env$upper = x
#}
#
#setLower = function(x,env){
#	env$lower = x
#}

#old function call commented out until memory issue is resolved
#e.divisive = function(X,sig.lvl=.05,R=199,eps=1e-3,half=1000,k=NULL,min.size=30,alpha=1,onlyOne=FALSE){
e.divisive = function(X,sig.lvl=.05,R=199,k=NULL,min.size=30,alpha=1){
  if(R < 0 && is.null(k))
  	stop("R must be a nonnegative integer.")
  if((sig.lvl <= 0 || sig.lvl >= 1) && is.null(k))
  	stop("sig.lvl must be a positive real number between 0 and 1.")
  if(min.size < 2)
  	stop("min.size must be an integer greater than 1.")
  if(alpha > 2 || alpha <= 0)
  	stop("alpha must be in (0,2].")
  n = nrow(X)#Length of time series
  energy = new.env(parent = emptyenv())
 
  if(is.null(k)){
	  k = n
	  #v = getBounds(R,sig.lvl,eps*(1:(R+1))/(1:(R+1)+half))
	  #setUpper(v$u,energy)
	  #setLower(v$l,energy)
	  #assign('upper',v$u,envir=energy) #energy$upper = v$u
	  #assign('lower',v$l,envir=energy) #energy$lower = v$l
  }
  else
	  R = 0 #no need to perform the permutation test

  ret = pvals = permutations = NULL
  changes = c(1,n+1)#list of change-points
  ret$k.hat = 1
  
  setDim(n,2,energy) #matrix used to avoid recalculating statistic for clusters
  D = as.matrix(dist(X))**alpha #distance matrix
  con = NULL 
  
  while(k>0){
    tmp = e.split(changes,D,min.size,FALSE,energy) #find location for change point
    i = tmp[[1]]; j = tmp[[2]]; Estat=tmp[[4]]; tmp = tmp[[3]]
    con = tail(tmp,1) #newest change-point
    if(con == -1)
      break #not able to meet minimum size constraint  
    result = sig.test(D,R,changes,min.size,Estat,env=energy) #run permutation test
    pval = result[1] #approximate p-value
    permutations = c(permutations,result[2]) #number of permutations performed
    pvals = c(pvals,pval) #list of computed p-values
    if(pval > sig.lvl)
      break #change point not significant
    changes = tmp #update set of change points
    ret$k.hat = ret$k.hat+1 #update number of clusters
    k = k-1
  }
  
  tmp = sort(changes)
  ret$order.found = changes
  ret$estimates = tmp
  ret$considered.last = con
  ret$p.values = pvals
  ret$permutations = permutations
  ret$cluster = rep(1:length(diff(tmp)),diff(tmp))
  return(ret)
}

e.split = function(changes,D,min.size,for.sim=FALSE, env=emptyenv()){ 	
 splits = sort(changes) #sort the set of current change points
 best = c(-1,-Inf)
 ii = jj = -1
 if(for.sim){ #If procedure is being used for significance test
	 for(i in 2:length(splits)){#iterate over intervals
		  tmp=splitPoint(splits[i-1],splits[i]-1,D,min.size)
		  if(tmp[2]>best[2]){ #tmp[2] is the "energy released" when the cluster was split
		    ii=splits[i-1]; jj=splits[i]-1
		    best=tmp #update the best point to split found so far
		  }
		}
		
		changes=c(changes,best[1]) #update the list of changepoints
		return(list('first'=ii,'second'=jj,'third'=changes,'fourth'=best[2]))
	}
	else{
		for(i in 2:length(splits)){ #iterate over intervals
			if(env$done_[splits[i-1],1])
				tmp = env$done_[splits[i-1],]
			else{
				tmp = splitPoint(splits[i-1],splits[i]-1,D,min.size)
				setVal(splits[i-1],1,tmp[1],env)
				setVal(splits[i-1],2,tmp[2],env)
			}
			if(tmp[2]>best[2]){ #tmp[2] is the "energy released" when the cluster was split
				ii = splits[i-1]
				jj = splits[i]-1
				best = tmp
			}
			
		}
		changes = c(changes,best[1]) #update the list of changepoints
		setVal(ii,1,0,env)#update matrix to account for newly proposed change point
		setVal(ii,2,0,env)#update matrix to account for newly proposed change point
		#energy$done_[ii,] = c(0,0) 
		return(list('first'=ii,'second'=jj,'third'=changes,'fourth'=best[2]))
	}
}

splitPoint = function(start,end,D,min.size){
	if(end-start+1 < 2*min.size) #interval too small to split
		return(c(-1,-Inf))
	D = D[start:end,start:end] #use needed part of distance matrix
	return(splitPointC(start,end,D,min.size))
}

sig.test = function(D,R,changes,min.size,obs,env=emptyenv()){
	if(R == 0) #no permutations so return a p-value of 0
		return(0)
	over = 0
	for(f in 1:R){
		D1=perm.cluster(D,changes) #permute within clusters
		tmp=e.split(changes,D1,min.size,TRUE)
		if(tmp[[4]] >= obs)
			over = over+1
		#if(over >= env$upper[f] || over <= env$lower)
		#	break
	}
		p.val = (1+over)/(f+1)
	return(c(p.val,f))
}

perm.cluster = function(D,points){
  points = sort(points)
  K = length(points)-1 #number of clusters
  for(i in 1:K ){ #shuffle within clusters by permuting matrix columns and rows
  	 u = sample(points[i]:(points[i+1]-1))
    D[points[i]:(points[i+1]-1),points[i]:(points[i+1]-1)] = D[u,u]
  }
  return(D)
}

