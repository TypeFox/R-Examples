#The only function that the user has to interact with. member is the initial membership
#array, Data is the observed data set, alpha is the exponent used on the euclidian 
#distance.
e.agglo = function(X,member=1:nrow(X),alpha=1,penalty = function(cps){0}){
	if(alpha<=0 || alpha>2)
		stop("The alpha argument must be in (0,2].")
	if(!is.function(penalty))
		stop("The penalty argument must be a function.")
	ret = process.data(member,X,alpha)
	N = ret$N
	for(K in N:(2*N-2)){
		#find which clusters optimize the GOF and then update the distnaces
		best = find.closest(K,ret)#find which clusters to merge
		ret$fit = c(ret$fit,best[3])#update GOF statistic
		ret = updateDistance(best[1],best[2],K,ret)#update information after merger
	}
	#penalize the GOF statistic
	cps = apply(ret$progression,1,function(x){x[!is.na(x)]})
	ret$fit = ret$fit + sapply(cps,penalty)
	#get the set of change points for the "best" clustering
	ret$estimates = which.max(ret$fit)
	ret$estimates = sort(ret$progression[ret$estimates,])
	#remove change point N+1 if a cyclic merger was performed
	if(ret$estimates[1] != 1)
		ret$estimates = rev(rev(ret$estimates)[-1])
	#create final membership vector
	tmp = ret$estimates
	if(tmp[1] == 1)
		ret$cluster = rep(1:length(diff(tmp)),diff(tmp))
	else{
		tmp = c(1,tmp)
		ret$cluster = rep(1:length(diff(tmp)),diff(tmp))
		k = nrow(X) - length(ret$cluster)
		ret$cluster = c(ret$cluster,rep(1,k))
	}
	#remove unnecessary output info
	ret$N = ret$left = ret$open = ret$D = ret$right = ret$lm = ret$sizes = NULL
	return(ret)
}

#Initialize all the necessary components of the list.
process.data = function(member,X,alpha){
	ret = NULL #the list with the necessary information
	u = unique(member)
	N = ret$N = length(u)#number of clusters
	fit = 0
	for(i in 1:N)#relabel clusters
		member[member==u[i]] = i
	if(is.unsorted(member))#Check that segments consist only of adjacent observations
		stop("Segments must be contiguous.")
	ret$sizes = numeric(2*N)
	ret$right = ret$left = array(0,dim=2*N-1)#left and right neighbor of a cluster
	ret$open = array(TRUE,dim=2*N-1)#TRUE means that a cluster has not been merged
	for(i in 1:N)#calculate initial cluster sizes
		ret$sizes[i] = sum(member==i)
	#set up left and right neighbors
	for(i in 2:(N-1)){
		ret$left[i] = i-1
		ret$right[i] = i+1
	}
	#special case for clusters 1 and N to allow for cyclic merging
	ret$left[N] = N-1
	ret$right[N] = 1
	ret$left[1] = N
	ret$right[1] = 2
	#matrix to say which clusters were merged at eaach step
	ret$merged = matrix(NA,nrow=N-1,ncol=2)
	#array of within distances
	within = numeric(N)
	for(i in 1:N)
		within[i] = getWithin(alpha,matrix(X[member==i,],ret$sizes[i],ncol(X)))
	#make distance matrix
	ret$D = matrix(Inf,nrow=2*N,ncol=2*N)
	for(i in 1:N)
		for(j in i:N)
			if(j != i)
				ret$D[i,j] = ret$D[j,i] = getBetween(alpha,matrix(X[member==i,],ret$sizes[i],ncol(X)),matrix(X[member==j,],ret$sizes[j],ncol(X))) - within[i] - within[j]
	diag(ret$D) = 0
	#set initial GOF value
	for(i in 1:N)
		fit = fit + ret$D[i,ret$left[i]] + ret$D[i,ret$right[i]]
	ret$fit = fit
	#create matrix for change point progression
	ret$progression = matrix(NA,nrow=N,ncol=N+1)
	ret$progression[1,1]=1
	for(i in 2:(N+1))
		ret$progression[1,i] = ret$progression[1,i-1] + ret$sizes[i-1]
	#vector to specify the starting point of a cluster
	ret$lm = numeric(2*N-1)
	ret$lm[1:N] = 1:N
	return(ret)
}

#Determine which clusters will be merged. Returns a vector of length 3. The first 
#element of the vector is the left cluster. The third element of the vector is the 
#newest GOF value.
find.closest = function(K,ret){
	N = ret$N
	best = -Inf
	triple = c(0,0,0)
	#iterate to see how the GOF value changes
	for(i in 1:K){
		if(ret$open[i]){
			x = gof.update(i,ret)#get updated gof value
			if(x > best){#better value found so update
				best = x
				triple = c(i,ret$right[i],x)
			}
		}
	}
	if(!triple[1] & !triple[2] & K!=(2*N-2))
		stop("There was a problem finding which clusters to merge.")
	return(triple)
}

#Function to calculate the new GOF value. The i argument is assumed to be the left 
#cluster.
gof.update = function(i,ret){
	fit = tail(ret$fit,1)
	j = ret$right[i]
	#get new left and right clusters
	rr = ret$right[j]
	ll = ret$left[i]
	#remove unneeded values in the GOF
	fit = fit - 2*(ret$D[i,j] + ret$D[i,ll] + ret$D[j,rr])
	#get cluster sizes
	n1 = ret$sizes[i]
	n2 = ret$sizes[j]
	#add distance to new left cluster
	n3 = ret$sizes[ll]
	k = ((n1+n3)*ret$D[i,ll] + (n2+n3)*ret$D[j,ll] - n3*ret$D[i,j])/(n1+n2+n3)
	fit = fit + 2*k
	#add distance to new right cluster
	n3 = ret$sizes[rr]
	k = ((n1+n3)*ret$D[i,rr] + (n2+n3)*ret$D[j,rr] - n3*ret$D[i,j])/(n1+n2+n3)
	fit = fit + 2*k
	return(fit)
}

#Function to update the distance from the new cluster to the other clusters. 
#i is assumed to be the left cluster. Also #updates any other necessary 
#information in the list (ret)
updateDistance = function(i,j,K,ret){
	#say which clusters were merged
	if(i <= ret$N)
		ret$merged[K-ret$N+1,1] = -i
	else
		ret$merged[K-ret$N+1,1] = i-ret$N
	if(j <= ret$N)
		ret$merged[K-ret$N+1,2] = -j
	else
		ret$merged[K-ret$N+1,2] = j-ret$N
	#update left and right neighobors
	ll = ret$left[i]
	rr = ret$right[j]
	ret$left[K+1] = ll
	ret$right[K+1] = rr
	ret$right[ll] = ret$left[rr] = K+1
	#update information on which clusters have been merged
	ret$open[i] = ret$open[j] = FALSE
	#assign size to newly created cluster
	n1 = ret$sizes[i]
	n2 = ret$sizes[j]
	ret$sizes[K+1] = n1+n2
	#update set of change points
	ret$progression[K-ret$N+2,] = ret$progression[K-ret$N+1,]
	ret$progression[K-ret$N+2,ret$lm[j]] = NA
	ret$lm[K+1] = ret$lm[i]
	#update distnaces
	for(k in 1:K){
		if(ret$open[k]){
			n3 = ret$sizes[k]
			n = n1+n2+n3
			ret$D[K+1,k] = ret$D[k,K+1] = ((n-n2)*ret$D[i,k]+(n-n1)*ret$D[j,k]-n3*ret$D[i,j])/n
		}
	}
	return(ret)
}
