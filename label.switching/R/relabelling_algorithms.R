# ECR algorithm:
# zpivot: the pivot
# z: mcmc output of allocation variables with dimension equal to mxn (m denotes MCMC iterations)
ecr<-function(zpivot,z,K){
if(K<max(z)){stop("K should be at least equal to max(z)")}
if(dim(z)[2]!=length(zpivot)){if(K<max(z)){stop("zpivot has not the same length with simulated z's")}}
n<-length(zpivot)
m<-dim(z)[1]
k<-K
# ECR algorithm
# By considering the problem of maxizing the S similarity between z[,iter] and zpivot
# as a special case of the assignment problem, it minimizes the corresponding cost matrix 
# between z[,iter] and zpivot
#library(lpSolve);				# contains the routine lp.assign for the solution of the assignment problem
st <- 1:k;							# the set {1,...,k}
perm <- array(data = NA, dim=c(m,K));	# the permutation of st that will be used to reorder the parameter values at each iteration
cost.matrix <- matrix(numeric(k*k),nrow=k,ncol=k); # (k times k) cost matrix of the assignment problem
s <- 1:n
#ptm<-proc.time()
for(iter in 1:m){
	alloc <- z[iter,]
	for(i in 1:k){
		ind <- which(alloc==i)
		so <- zpivot[ind];	# finding the indices of zpivot that correspond to z[,iter]==i
		l <- length(ind);	# the length of the vector above
		for(j in 1:k) cost.matrix[i,j] <- l-length(so[so==j])		# the cost of assigning to index i the permutation j
		}
	matr <- lp.assign(cost.matrix)$solution ;	# solution: matr[j,i] = 1 <=> index i assigned to index j
	for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]}	# the optimal permutation for the current iteration
	}
	#time<-proc.time() - ptm
#	print(paste("ECR algorithm time:", as.numeric(time[3]),"seconds."))
	results<-list(perm)
	names(results)<-c("permutations")
	return(results)
}
#################################################################################
#################################################################################
ecr.iterative.1<-function (z, K, opt_init,threshold,maxiter) {
    if (K < max(z)) {
        stop("K should be at least equal to max(z)")
    }
    m <- dim(z)[1]
    k <- K
    st <- 1:k
    perm <- array(data = NA, dim = c(m, K))
    cost.matrix <- matrix(numeric(k * k), nrow = k, ncol = k)
    n <- dim(z)[2]
    s <- 1:n
    if (missing(opt_init)) {
        for (j in 1:K) {
            perm[, j] <- j
        }
    }
    else {
        perm <- opt_init
    }
    permz <- z
    for(iter in 1:m){
	permz[iter, ] <- perm[iter, ][z[iter, ]]
    }	
    zpivot <- numeric(n)
    criterion <- 99
#    threshold <- 10^(-6)
    if(missing(threshold)){threshold <- 10^(-6)}
    if(missing(maxiter)){maxiter <- 100}
    t <- 1
    zpivot <- apply(permz, 2, function(y) order(table(y), decreasing = T)[1])
    cf <- 0
    for (iter in 1:m) {
        alloc <- z[iter, ]
        for (i in 1:k) {
            so <- zpivot[s[alloc == i]]
            l <- length(alloc[alloc == i])
            for (j in 1:k) cost.matrix[i, j] <- l - length(so[so == 
                j])
        }
        matr <- lp.assign(cost.matrix)$solution
        for (i in 1:k) {
            perm[iter, i] <- st[matr[, i] > 0]
        }
        cf <- cf + sum(cost.matrix * matr)
        permz[iter, ] <- perm[iter, ][z[iter, ]]
    }
    previous <- cf
    #print(paste("t = ",t,"cost function = ", cf))
    while ((criterion > threshold)&&(t<maxiter)) {
	previousperm<-perm
  #   while (t < 4) {
        t <- t + 1
        zpivot <- apply(permz, 2, function(y) order(table(y), 
            decreasing = T)[1])
        cf <- 0
        for (iter in 1:m) {
            alloc <- z[iter, ]
            for (i in 1:k) {
                so <- zpivot[s[alloc == i]]
                l <- length(alloc[alloc == i])
                for (j in 1:k) cost.matrix[i, j] <- l - length(so[so == 
                  j])
            }
            matr <- lp.assign(cost.matrix)$solution
            for (i in 1:k) {
                perm[iter, i] <- st[matr[, i] > 0]
            }
            cf <- cf + sum(cost.matrix * matr)
            permz[iter, ] <- perm[iter, ][z[iter, ]]
        }
        current <- cf
        criterion <- previous - current
        previous <- cf
	if (criterion<0){
		perm<-previousperm
	}
        #print(paste("t = ",t,"cost function = ", cf))
    }

    status <- paste("Converged (",t," iterations)",sep="")
    if(criterion>threshold){status <- "Max iterations exceeded"}

    results <- list(perm, t,status)
    names(results) <- c("permutations", "iterations","status")
    return(results)
}









###################################################################################
#######################################################################################
#######################################################################################
#######################################################################################





ecr.iterative.2<-function(z,K,p,threshold,maxiter){
if(K<max(z)){stop("K should be at least equal to max(z)")}
n<-dim(z)[2]
m<-dim(z)[1]
k<-K
#library(lpSolve);				# contains the routine lp.assign for the solution of the assignment problem
st <- 1:k;							# the set {1,...,k}
perm <- array(data = NA, dim=c(m,K));	# the permutation of st that will be used to reorder the parameter values at each iteration
cost.matrix <- matrix(numeric(k*k),nrow=k,ncol=k); # (k times k) cost matrix of the assignment problem
s <- 1:n
#ptm<-proc.time()

#initial permutations: identity
for (j in 1:K){
perm[,j]<-j
}


zpivot<-numeric(n)

criterion<-99
#threshold<-10**(-6)
if(missing(threshold)){threshold <- 10^(-6)}
if(missing(maxiter)){maxiter <- 100}
t<-1
#estimating pivot

q<-array(data = 0, dim =c(n,k))
for (j in 1:k){
	for(iter in 1:m){q[,j]<-q[,j] + p[iter,,perm[iter,j]]}
}
q<-q/m

for(i in 1:n){zpivot[i]<-order(q[i,])[K]}
cf<-0
for(iter in 1:m){
	alloc <- z[iter,]
	for(i in 1:k){
		so <- zpivot[s[alloc==i]];	# finding the indices of zpivot that correspond to z[,iter]==i
		l <- length(alloc[alloc==i]);	# the length of the vector above
		for(j in 1:k)cost.matrix[i,j] <- l-length(so[so==j])
	}
	matr <- lp.assign(cost.matrix)$solution ;	# solution: matr[j,i] = 1 <=> index i assigned to index j
	for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]}	# the optimal permutation for the current iteration
	cf<-cf + sum(cost.matrix*matr)
}
previous<-cf
#print(paste("iteration 1, cost function =",cf))

while((criterion>threshold)&&(t<maxiter)){
	t<-t + 1
	#estimating pivot
	q<-array(data = 0, dim =c(n,k))
	for (j in 1:k){
		for(iter in 1:m){q[,j]<-q[,j] + p[iter,,perm[iter,j]]}
	}
	q<-q/m

	for(i in 1:n){zpivot[i]<-order(q[i,])[K]}
	cf<-0
	for(iter in 1:m){
		alloc <- z[iter,]
		for(i in 1:k){
			so <- zpivot[s[alloc==i]];	# finding the indices of zpivot that correspond to z[,iter]==i
			l <- length(alloc[alloc==i]);	# the length of the vector above
			for(j in 1:k)cost.matrix[i,j] <- l-length(so[so==j])
		}
		matr <- lp.assign(cost.matrix)$solution ;	# solution: matr[j,i] = 1 <=> index i assigned to index j
		for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]}	# the optimal permutation for the current iteration
		cf<-cf + sum(cost.matrix*matr)
	}
	current<-cf
	criterion<-abs(previous - current)
	previous<-cf
	#print(paste("iteration", t, "criterion =",cf))


}
	#time<-proc.time() - ptm
	#print(paste("Iterative ECR algorithm 2 converged at", t,"iterations"))
	#print(paste("Iterative ECR algorithm 2 time:", as.numeric(time[3]),"seconds."))
	status <- paste("Converged (",t," iterations)",sep="")
	if(criterion>threshold){status <- "Max iterations exceeded"}

	results<-list(perm,t,status)
	names(results)<-c("permutations","iterations","status")
	return(results)
}












###################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
# p: mxnxK matrix of allocation proababilities

stephens<-function(p,threshold,maxiter){
n<-dim(p)[2]
m<-dim(p)[1]
k<-dim(p)[3]
burnin<-0
K<-k
#library(lpSolve);				# contains the routine lp.assign for the solution of the assignment problem
st <- 1:k;							# the set {1,...,k}
perm <- c(numeric(k));	# the permutation of st that will be used to reorder the parameter values at each iteration
cost.matrix <- matrix(numeric(k*k),nrow=k,ncol=k); # (k times k) cost matrix of the assignment problem
s <- 1:n
###

up.threshold<-1-10**(-6)
down.threshold<-10**(-6)
for(k in 1:K){
for(i in 1:n){
up.index<-which(p[,i,k]>up.threshold)
down.index<-which(p[,i,k]<down.threshold)
if(length(up.index)>0){p[up.index,i,k]<-rep(up.threshold,length(up.index))}
if(length(down.index)>0){p[down.index,i,k]<-rep(down.threshold,length(down.index))}
}
}
for(iter in 1:m){
p[iter,,]<-p[iter,,]/rowSums(p[iter,,])
}
##
#print(paste("finished smoothing"))
#step 0.
perm <- array(data = NA, dim = c(m,k))
for(j in 1:k){perm[,j]<-j}
q<-array(data = 0, dim =c(n,k))
previous<- -99
criterion<-99
#threshold<-10**(-6)
if(missing(threshold)){threshold <- 10^(-6)}
if(missing(maxiter)){maxiter <- 100}
#maxiter<-99
#ptm<-proc.time()
# t denotes stephens algorithm iterations
t<-0
while((criterion>threshold)&&(t<maxiter)){
	t<-t+1
#compute q matrix
	q<-array(data = 0, dim =c(n,k))
	for (j in 1:k){
		for(iter in 1:m){q[,j]<-q[,j] + p[iter,,perm[iter,j]]}
	}
	q<-q/m
	
	for(iter in 1:m){
		for(j in 1:k){		
			temp<-p[iter,,]*(log(p[iter,,]) - log(q[,j]))
			cost.matrix[j,]<-colSums(temp)
		}
		matr <- lp.assign(cost.matrix)$solution
		for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]
		}
		perm[iter,]<-order(perm[iter,])
	}
	current<-cost.function<-sum(cost.matrix)
	criterion<-abs(previous - current)
	previous<-current
	#if(t>1){print(paste("iteration ", t,", criterion =",criterion))}
#this ends t
#
}
	#if(t<maxiter){print(paste("Stephens' algorithm converged at", t,"iterations"))}else{
	#	print(paste("Stephens' algorithm didn't converged at", maxiter,"iterations"))
	#}
	status <- paste("Converged (",t," iterations)",sep="")
	if(criterion>threshold){status <- "Max iterations exceeded"}
	#time<-proc.time() - ptm
	#print(paste("Stephens' algorithm time:", as.numeric(time[3]),"seconds."))
	results<-list(perm,t,status)
	names(results)<-c("permutations","iterations","status")
	return(results)
}
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


sjw<-function(mcmc.pars,z,complete,x,init,threshold,maxiter){
K<-dim(mcmc.pars)[2]
J<-dim(mcmc.pars)[3]
m<-dim(mcmc.pars)[1]
ptm <- proc.time()

#initial estimate
# if init = 0 algorithm starts from the overall mean.
if(init <1){
#print(paste("Initialization using the overall mean"))
theta.hat<-array(data = NA, dim =c(K,J))
for(k in 1:K){
for(j in 1:J){
theta.hat[k,j]<-mean(mcmc.pars[,k,j])
}
}
}else{
init = floor(init + 0.5)
#print(paste("Initialization using MCMC iteration: ",init))
theta.hat<-mcmc.pars[init,,]}
#print(paste(theta.hat))
previous<-theta.hat
#library(combinat)
permutations<-permn(1:K)
kfactorial<-gamma(K+1)
#permutation probabilities
pp<-array(data = NA, dim =c(m,kfactorial))
#t denotes EM iterations
#maxiter<-100
t<-0
#threshold<-10**(-4)
if(missing(threshold)){threshold <- 10^(-6)}
if(missing(maxiter)){maxiter <- 100}
criterion <- 99
while((criterion>threshold)&&t<maxiter){
t<-t+1
temp<-numeric(kfactorial)
#E-step: compute the permutation probabilities
for(iter in 1:m){
for(k in 1:kfactorial){
perm<-permutations[[k]]
temp[k]<-complete(x,perm[z[iter,]],theta.hat)
}

for(k in 1:(kfactorial-1)){
pp[iter,k]<-1/sum(exp(temp-temp[k]))
}
pp[iter,kfactorial]<- 1 - sum(pp[iter,1:(kfactorial-1)])
if(is.na(max(pp[iter,]))==TRUE){pp[iter,]=rep(0,kfactorial);ind<-order(temp)[kfactorial];pp[iter,ind]<-1}
index<-which(pp[iter,]<0)
if(length(index)>0){pp[iter,index]=rep(0,length(index))}

}

#M-step: update parameters estimate
theta.hat<-array(data = 0, dim =c(K,J))
for(iter in 1:m){
for(k in 1:kfactorial){
for(j in 1:J){
perm<-permutations[[k]]
theta.hat[,j]<-theta.hat[,j] + pp[iter,k]*mcmc.pars[iter,perm,j]
}
}
}

theta.hat<-theta.hat/(m)
current<-theta.hat
criterion<-max(abs(previous - current))
previous<-theta.hat
#cat(paste("iter = ",t, ", criterion = ",criterion),'\n')
#print(paste(theta.hat))


#this ends t

}

status <- paste("Converged (",t," iterations)",sep="")
if(criterion>threshold){status <- "Max iterations exceeded"}
perm <- array(data = NA, dim=c(m,K));
#sample of permutations per iteration
for (i in 1:m){
k<-sample(1:kfactorial,1,prob = pp[i,])
perm[i,]<-order(permutations[[k]])
}

#time<-proc.time() - ptm
#print(paste("SJW algorithm time:", as.numeric(time[3]),"seconds."))

results<-list(perm,t,status)
names(results)<-c("permutations","iterations","status")
return(results)

}



##########################################################################################
##########################################################################################
#		PRA
pra<-function(mcmc.pars,pivot){
K<-dim(mcmc.pars)[2]
J<-dim(mcmc.pars)[3]
m<-dim(mcmc.pars)[1]
#if(K>8){cat(paste("WARNING: PRA is not suggested for", K,"components"),'\n')}
l1<-dim(pivot)
if((l1[1]!=K)||(l1[2]!=J)){
stop("Pivot and MCMC samples should have the same length")
}
perms <- array(data = NA, dim = c(m,K))
#library(combinat)
#ptm <- proc.time()
permutations<-permn(1:K)
kfactorial<-gamma(K+1)
for(iter in 1:m){
	t<-1
	perm<-permutations[t][[1]]
	inner.product<-sum(pivot*mcmc.pars[iter,perm,])
	max.val<-inner.product
	perms[iter,]<-perm
	for(t in 2:kfactorial){
		perm<-permutations[t][[1]]
		inner.product<-sum(pivot*mcmc.pars[iter,perm,])
		if(inner.product>max.val){perms[iter,]<-perm;max.val<-inner.product}
	}
}
#	time<-proc.time() - ptm
#	print(paste("PRA time:", as.numeric(time[3]),"seconds."))

	results<-list(perms)
	names(results)<-c("permutations")
	return(results)
}


############################################################################
############################################################################
############################################################################
###########################################################################



#		aic: artificial identifiability constraint
# by default the ordering constraint is imposed to the first parameter
# constrain = 1,...,J controls the position of the constraint
aic <- function(mcmc.pars,constraint){
K<-dim(mcmc.pars)[2]
J<-dim(mcmc.pars)[3]
m<-dim(mcmc.pars)[1]
if((constraint %in% 1:J)==FALSE){
	stop(cat(paste("constraint should be integer between 1 and",J),"\n"))	
}
perms <- array(data = NA, dim = c(m,K))
for(iter in 1:m){
	perms[iter,] <- order(mcmc.pars[iter,,constraint])
}
	results<-list(perms)
	names(results)<-c("permutations")
	return(results)
}

##########################################################################################
##########################################################################################


#data based algorithm (rodriguez and walker, jcgs)
dataBased <- function(x,K,z){
	# x: data
	# K: number of components
	if(is.null(dim(x))==TRUE)x<-array(x,dim=c(length(x),1))
	p <- dim(x)[2]
	n <- dim(x)[1]
	m <- dim(z)[1]
	if (dim(z)[2] != n) {
		stop("sample size has not the same length with simulated z's")
	}
	#step 1. data-based relabelling: finding estimates ("algorithm 5")
	cluster_mean <- array(data = 0, dim = c(K,p))
	cluster_sd <- array(data = 0, dim = c(K,p))
	cluster_n <- numeric(K)
	#initialize
	n_mean <- n_sd <- rep(1,K)
	for(j in 1:p){
		data_r <- diff(range(x[,j]))
		data_min <- min(x[,j])
		for(k in 1:K){
			cluster_mean[k,j] <- data_min + data_r*j/(k+1)			
			cluster_sd[k,j] <- data_r/K
		}		
	}
	#update estimates
	perm <- numeric(K)
	cost.matrix <- matrix(numeric(K * K), nrow = K, ncol = K)
	st <- 1:K
	s <- 1:n
	sample_mean <- sample_sd <- array(data = 0,dim=c(K,p))
	for(iter in 1:m){
		alloc <- z[iter, ]
		for (i in 1:K) {
			#so <- zpivot[s[alloc == i]]
			ind <- which(alloc == i)
			cluster_n[i] <- length(ind)
			for (j in 1:k){
				cost.matrix[i, j] <- 0
				for(g in 1:p){
					cost.matrix[i, j] <- cost.matrix[i, j] + sum(((x[ind,g] - cluster_mean[j,g])/cluster_sd[j,g])**2)
				}
			for(g in 1:p){
				sample_mean[i,g] <- sum(x[ind,g])/cluster_n[i]
				if(cluster_n[i]>1)sample_sd[i,g] <- sqrt(sum((x[ind,g]-sample_mean[i,g])**2)/(cluster_n[i]-1))
			} 

			}
		}
		matr <- lp.assign(cost.matrix)$solution
		for (i in 1:k){
			perm[i] <- st[matr[, i] > 0]
			permSum <- cluster_n[perm[i]]
			if(permSum>0){
				cluster_mean[i,] <- ((n_mean[i]-1)*cluster_mean[i,]+sample_mean[perm[i],])/n_mean[i]
				n_mean[i] <- n_mean[i] + 1
				if(permSum>1){
					cluster_sd[i,] <- ((n_sd[i]-1)*cluster_sd[i,]+sample_sd[perm[i],])/n_sd[i]
					n_sd[i] <- n_sd[i] + 1					
				}
			}
		}
	}

	#step 2. data-based relabelling: estimates strategy ("algorithm 6")
	perm <- array(data = NA, dim = c(m, K))
	for(iter in 1:m){
		alloc <- z[iter, ]
		for (i in 1:K) {
			ind <- which(alloc == i)
			cluster_n[i] <- length(ind)
			for (j in 1:k){
				cost.matrix[i, j] <- 0
				for(g in 1:p){
					cost.matrix[i, j] <- cost.matrix[i, j] + sum(((x[ind,g] - cluster_mean[j,g])/cluster_sd[j,g])**2)
				} 
			}
		}
		matr <- lp.assign(cost.matrix)$solution
		for (i in 1:k) {
            		perm[iter, i] <- st[matr[, i] > 0]
        	}
	}

	results <- list(perm)
	names(results) <- c("permutations")
	return(results)


}
###########################################################


############################################################################


permute.mcmc<-function(mcmc,permutations){
m<-dim(permutations)[1]
K<-dim(permutations)[2]
J<-dim(mcmc)[3]
mcmc.permuted<-mcmc
	for(iter in 1:m){
		for(j in 1:J){
			mcmc.permuted[iter,,j]<-mcmc[iter,permutations[iter,],j]
		}	
	}
	results<-list(mcmc.permuted)
	names(results)<-c("output")
	return(results)
}

################################################################################################
################################################################################################
################################################################################################



label.switching<-function(method, zpivot, z, K, prapivot, p, complete, mcmc, sjwinit, data, constraint, groundTruth,thrECR,thrSTE,thrSJW,maxECR,maxSTE,maxSJW,userPerm){
    runStart <- proc.time()
    cat('\n')
    L <- length(method)
    for (l in 1:L) {
        if ((method[l] %in% c("ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", 
            "STEPHENS", "SJW", "PRA","DATA-BASED","AIC","USER-PERM")) == FALSE) {
            stop(cat(paste("method:", method[l], "is not recognised"),'\n'))
        }
    }


    if (missing(z) == FALSE) {
        m <- dim(z)[1]
	nCheck <- dim(z)[2]
	kMinCheck <- max(z)
	kCheck <- kMinCheck
    }
    else {
        if (missing(mcmc) == FALSE) {
            m <- dim(mcmc)[1]
	    kCheck <- dim(mcmc)[2] 
	    kMinCheck <- kCheck
        }
        else {
            if (missing(p) == FALSE) {
                m <- dim(p)[1]
		kCheck <- dim(p)[3]
		nCheck <- dim(p)[2]
		kMinCheck <- kCheck
            }
            else {
                stop(cat(paste("    [ERROR]: At least one of z, mcmc, or p should be provided"),'\n'))
            }
        }
    }
if (missing(p) == FALSE){nCheck <- dim(p)[2]}
if (missing(p) == FALSE){kCheck <- dim(p)[3]}
if (missing(mcmc) == FALSE){kCheck <- dim(mcmc)[2]}


if (missing(data) == FALSE){
	if(is.null(dim(data))==TRUE){nCheck = length(data)}else{
		nCheck <- dim(data)[1]
	}
}
if(missing(K)==TRUE){K = kCheck;
cat(paste("    [WARNING]: K is not provided. According to input it is assumed that K = ",K,".",sep=""),'\n')
}
if(max(z) < K){
cat(paste("    [WARNING]: max sampled latent allocation = ",max(z)," < ", "K = ",K,". This indicates that the MCMC sampler has not converged and/or the presence of redundant components.",sep=""),'\n')
}


## checking input consistency
if(missing(z)==FALSE){
	if(m != dim(z)[1]){stop(cat(paste("    [ERROR]: MCMC iterations are not equal to dim(z)[1]."),'\n'))}
}
if(missing(mcmc)==FALSE){
	if(m != dim(mcmc)[1]){stop(cat(paste("    [ERROR]: MCMC iterations are not equal to dim(mcmc)[1]."),'\n'))}
}
if(missing(p)==FALSE){
	if(m != dim(p)[1]){stop(cat(paste("    [ERROR]: MCMC iterations are not equal to dim(p)[1]."),'\n'))}
}
if(missing(mcmc)==FALSE){
	if(missing(kCheck)){kCheck = dim(mcmc)[2]}
	if(kCheck != dim(mcmc)[2]){stop(cat(paste("    [ERROR]: K is not equal to dim(mcmc)[2]."),'\n'))}
}
if(missing(p)==FALSE){
	if(missing(kCheck)){dim(p)[3]}
	if(kCheck != dim(p)[3]){stop(cat(paste("    [ERROR]: K is not equal to dim(p)[3]."),'\n'))}
}
if(missing(K)==FALSE){
	if(kCheck != K){stop(cat(paste("    [ERROR]: Number of components is not consistent with the input."),'\n'))}
	if(kMinCheck > K){stop(cat(paste("    [ERROR]: Number of components should be at least equal to ",kMinCheck,",", sep=""),'\n'))}
}

if(missing(z)==FALSE){
	if(nCheck != dim(z)[2]){stop(cat(paste("    [ERROR]: Number of observations is not equal to dim(z)[2]."),'\n'))}
}
if(missing(p)==FALSE){
	if(nCheck != dim(p)[2]){stop(cat(paste("    [ERROR]: Number of observations is not equal to dim(p)[2]."),'\n'))}
}
if(missing(zpivot)==FALSE){
	if(is.null(dim(zpivot))==TRUE){
		if(nCheck != length(zpivot)){stop(cat(paste("    [ERROR]: Number of observations is not equal to length(zpivot)."),'\n'))}
	}else{
		if(nCheck != dim(zpivot)[2]){stop(cat(paste("    [ERROR]: Number of observations is not equal to dim(zpivot)[2]."),'\n'))}
	}
}
if(missing(prapivot)==FALSE){
	if(kCheck != dim(prapivot)[1]){stop(cat(paste("    [ERROR]: K is not equal to dim(prapivot)[1]."),'\n'))}
}
if(missing(prapivot)==FALSE){
	if(dim(mcmc)[3] != dim(prapivot)[2]){stop(cat(paste("    [ERROR]: J is not equal to dim(prapivot)[2]."),'\n'))}
}
if(missing(data)==FALSE){
	if(is.null(dim(data))==TRUE){nX <- length(data)}else{
		nX <- dim(data)[1]
	}
	if(nCheck != nX){stop(cat(paste("    [ERROR]: data length is not compatible with the input."),'\n'))}
}
if(missing(groundTruth)==FALSE){
	if(length(groundTruth) != nCheck){stop(cat(paste("    [ERROR]: length(groundTruth) is not equal to number of observations."),'\n'))}
	if( all(groundTruth == floor(groundTruth)) == FALSE ){stop(cat(paste("    [ERROR]: non-integer groundTruth entries are not allowed."),'\n'))}
}
if(("ECR"%in%method)==TRUE){
	if(missing(z)==TRUE)stop(cat(paste("    [ERROR]: z is required for ECR."),'\n'))
	if(missing(zpivot)==TRUE)stop(cat(paste("    [ERROR]: zpivot is required for ECR."),'\n'))
}
if(("ECR-ITERATIVE-1"%in%method)==TRUE){
	if(missing(z)==TRUE)stop(cat(paste("    [ERROR]: z is required for ECR-ITERATIVE-1."),'\n'))
}
if(("ECR-ITERATIVE-2"%in%method)==TRUE){
	if(missing(z)==TRUE)stop(cat(paste("    [ERROR]: z is required for ECR-ITERATIVE-2."),'\n'))
}
if(("ECR-ITERATIVE-2"%in%method)==TRUE){
	if(missing(p)==TRUE)stop(cat(paste("    [ERROR]: p is required for ECR-ITERATIVE-2."),'\n'))
}
if(("STEPHENS"%in%method)==TRUE){
	if(missing(p)==TRUE)stop(cat(paste("    [ERROR]: p is required for STEPHENS."),'\n'))
}
if(("SJW"%in%method)==TRUE){
	if(missing(data)==TRUE)stop(cat(paste("    [ERROR]: data is required for SJW."),'\n'))
}
if(("SJW"%in%method)==TRUE){
	if(missing(complete)==TRUE)stop(cat(paste("    [ERROR]: complete is required for SJW."),'\n'))
}
if(("SJW"%in%method)==TRUE){
	if(missing(z)==TRUE)stop(cat(paste("    [ERROR]: z is required for SJW."),'\n'))
}
if(("SJW"%in%method)==TRUE){
	if(missing(mcmc)==TRUE)stop(cat(paste("    [ERROR]: mcmc is required for SJW."),'\n'))
}
if(("AIC"%in%method)==TRUE){
	if(missing(mcmc)==TRUE)stop(cat(paste("    [ERROR]: mcmc is required for AIC."),'\n'))
}
if(("DATA-BASED"%in%method)==TRUE){
	if(missing(z)==TRUE)stop(cat(paste("    [ERROR]: z is required for DATA-BASED."),'\n'))
}
if(("DATA-BASED"%in%method)==TRUE){
	if(missing(data)==TRUE)stop(cat(paste("    [ERROR]: data is required for DATA-BASED."),'\n'))
}

if(kCheck < 2){stop(cat(paste("    [ERROR]: K should be at least equal to 2."),'\n'))}
if(nCheck < 2){stop(cat(paste("    [ERROR]: n should be at least equal to 2."),'\n'))}
if(("PRA" %in% method)&&(kCheck>8)){cat(paste("    [WARNING]: PRA is not suggested for ", kCheck,"components"),'\n')}
if(("SJW" %in% method)&&(kCheck>8)){cat(paste("    [WARNING]: SJW is not suggested for ", kCheck,"components"),'\n')}
if(missing(thrECR)){thrECR <- 10^(-6)}
if(missing(thrSTE)){thrSTE <- 10^(-6)}
if(missing(thrSJW)){thrSJW <- 10^(-6)} 
if(missing(maxECR)){maxECR <- 100}
if(missing(maxSTE)){maxSTE <- 100}
if(missing(maxSJW)){maxSJW <- 100} 
minThreshold <- 1e-12
thrECR <- max(minThreshold,thrECR)
thrSTE <- max(minThreshold,thrSTE)
thrSJW <- max(minThreshold,thrSJW)
if(maxECR < 1){maxECR <- 100}
if(maxSTE < 1){maxSTE <- 100}
if(maxSJW < 1){maxSJW <- 100} 
    if(("USER-PERM" %in% method)==TRUE){
	if(missing(K)==TRUE){stop(cat(paste("   [ERROR]: K is not supplied."),'\n'))}
	#cat(paste("   [WARNING]: User-defined permutations supplied."),'\n')
	if(is.list(userPerm)==FALSE){
		temp <- vector('list',length=1)
		temp[[1]]<-userPerm
		userPerm <- temp
		temp <- 0
		userLength <- 1
	}else{
		userLength <- length(userPerm)
		if(is.null(names(userPerm))==TRUE){
			names(userPerm) <- paste("user",(1:userLength),sep="-")
		}
	}
	names(userPerm) <- paste("USER",(1:userLength),sep="-")
	cat(paste('    Checking user-supplied permutations for consistency...'))
	userStatus <- rep(1,userLength)
	for(j in 1:userLength){
		if(dim(z)[1] != dim(userPerm[[j]])[1]){userStatus[j]=0;cat('\n');cat(paste("    [ERROR]: number of MCMC samples should be equal to number of permutations."))}
		sss <- 1:K
		if(K != dim(userPerm[[j]])[2]){userStatus[j]=0;cat('\n');cat(paste("    [ERROR]: Number of components should be equal to permutation size."))}
		for(i in 1:dim(userPerm[[j]])[1]){
			myCheck <- table(match(unique(userPerm[[j]][i,]),sss))
			if(length(myCheck) != K){
				userStatus[j]=0;
				cat(paste('\n'))
				cat(paste('    problem in line ',i,' of supplied permutation set ',j,':',sep=""),'\n')
				cat(paste('   '),paste(userPerm[[j]][i,]),paste('is not a permutation of {1,...,',K,'}'),'\n')
				stop(cat(paste("    [ERROR]: user-defined input is not valid"),'\n'))
			}
		}
		#cat(paste('        permutation set ',j,':      OK',sep=""),'\n')
	}
	cat(paste('done.'),'\n')
    }


##############################
    fr <- 1
    dimname <- L
    if ((is.array(zpivot) == TRUE) && (("ECR" %in% method) == 
        TRUE)) {
        fr <- dim(zpivot)[1]
        dimname <- L + fr - 1
    }
    if ("ECR" %in% method) {
        ind <- which(method == "ECR")
        if (length(ind) > 1) {
            stop(paste("ECR appearing more than 1 times"))
        }
        nam <- numeric(dimname)
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind - 1)) {
                nam[i] <- method[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- paste("ECR", i, sep = "-")
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- method[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "ECR"
        }
    }
    else {
        nam <- method
    }
 #   cat(nam,'\n')
 #   cat(dimname,'\n')
    nams1 <- nam
##
    if((missing(constraint)==FALSE) && (length(constraint) == 1) && (constraint == "ALL")){J <- dim(mcmc)[3];constraint = 1:J}
    fr <- 1
    if ((missing(constraint)==FALSE) && (length(constraint)>1) && (("AIC" %in% method) == TRUE)) {
        fr <- length(constraint)
        dimname <- dimname + fr - 1
    }
    if ("AIC" %in% method) {
	
        ind <- which(nams1=="AIC")
        if (length(ind) > 1) {
            stop(paste("AIC appearing more than 1 times"))
        }
        #nam <- numeric(dimname)
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind  - 1)) {
                nam[i] <- nams1[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- paste("AIC", i, sep = "-")
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- nams1[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "AIC"
        }
    }

    #cat(nam,'\n')
##
    nams2 <- nam
    if ("USER-PERM" %in% method) {
        fr <- userLength
        dimname <- dimname + fr - 1
	
        ind <- which(nams2=="USER-PERM")
        if (length(ind) > 1) {
            stop(paste("USER-PERM appearing more than 1 times"))
        }
        #nam <- numeric(dimname)
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind  - 1)) {
                nam[i] <- nams2[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- names(userPerm)[i] #paste("USER", i, sep = "-")
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- nams2[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "USER"
        }
    }

    #cat(paste(nam),'\n')
   # stop(cat(paste("here")))




myPrettyPrint <- function(gap0,gap1,gap2,gap3,word1,word2,word3){
	nch1 = nchar(as.character(word1));
	nch2 = nchar(as.character(word2));
	nch3 = nchar(as.character(word3));
	cat(rep("",gap0),".",word1,paste(rep("",max(gap1 - nch1,1))),word2,rep("",max(1,gap2-nch2)),word3,rep("",max(1,gap3 - nch3)),".",'\n')
}


##




    permutations <- vector("list", length = dimname)
    names(permutations) <- nam
    timings <- numeric(dimname)
    names(timings) <- nam
    f <- 0
    fic <- 0
    fuser <- 0	
    t <- 1
    userIterator <- 0
    gap0 <- 4
    gap1 <- 30
    gap2 <- 20
    gap3 <- 30
    #cat(paste("    -------------------------------------------\n"))
    #cat(paste("    Method                           Time (sec)\n"))
    #cat()
    cat(paste("    ......................................................................................\n"))
    myPrettyPrint(gap0,gap1,gap2,gap3,"Method","Time (sec)","Status    ")
    cat(paste("    ......................................................................................\n"))
    for (l in 1:L) {
        #cat(paste("   ", method[l]))
        if (method[l] == "ECR") {
            if (is.array(zpivot) == TRUE) {
                while (f < dim(zpivot)[1]) {
                  f <- f + 1
                  #if (f == 1) {
                   # cat(paste(" (using pivot", f, "of", dim(zpivot)[1]), 
                   #   ")", sep = "")
                  #}
                  #else {
                  #  cat(paste("        (using pivot", f, "of", dim(zpivot)[1]), ")", sep = "")
                  #}
                  if (dim(zpivot)[2] != dim(z)[2]) {
                    stop(paste("length(zpivot) and number of columns of z are not equal"))
                  }
                  if (K < max(z)) {
                    stop(paste("K should be at least equal to", 
                      max(z)))
                  }
                  tpm <- proc.time()
                  permutations[[t]] <- ecr(zpivot[f, ], z, K)$permutations
                  time <- proc.time() - tpm
                  time <- round(as.numeric(time[3]),3) #as.numeric(time[3])
                  #cat(paste("        ", time, "\n"))
		  voutsas <- paste(method[l]," (pivot ", f, " of ", dim(zpivot)[1],")",sep="")
		  myPrettyPrint(gap0,gap1,gap2,gap3,voutsas,time,"OK")
                  timings[t] <- time
                  t <- t + 1
                }
            }
            else {
                if (missing(zpivot)) {
                  stop(paste("zpivot is missing"))
                }
                else {
                  if (length(zpivot) != dim(z)[2]) {
                    stop(paste("length(zpivot) and number of columns of z are not equal"))
                  }
                  if (K < max(z)) {
                    stop(paste("K should be at least equal to", 
                      max(z)))
                  }
                  tpm <- proc.time()
                  permutations[[t]] <- ecr(zpivot, z, K)$permutations
                  time <- proc.time() - tpm
                  time <- round(as.numeric(time[3]),3) #as.numeric(time[3])
                  #cat(paste("                             ", time, "\n"))
		  myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,"OK")
                  timings[t] <- time
                  t <- t + 1
                }
            }
        }
        if (method[l] == "ECR-ITERATIVE-1") {
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- ecr.iterative.1(z, K,threshold = thrECR,maxiter = maxECR)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]),3) # as.numeric(time[3])
            #cat(paste("                 ", time, "  (converged at", hold$iterations, "iterations)\n"))
	    myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "ECR-ITERATIVE-2") {
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(p)) {
                stop(paste("p is missing"))
            }
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- ecr.iterative.2(z, K, p,thrECR,maxECR)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]),3) #as.numeric(time[3])
            #cat(paste("                 ", time, "  (converged at", hold$iterations, "iterations)\n"))
	    myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "PRA") {
            tpm <- proc.time()
            permutations[[t]] <- pra(mcmc, prapivot)$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]),3) #as.numeric(time[3])
            #cat(paste("                             ", time, "\n"))
	    myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,"OK")
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "STEPHENS") {
            if (missing(p)) {
                stop(paste("p is missing"))
            }
            tpm <- proc.time()
            hold <- stephens(p,thrSTE,maxSTE)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]),3) #as.numeric(time[3])
            #cat(paste("                        ", time, "  (converged at", hold$iterations, "iterations)\n"))
	    myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "SJW") {
            if (missing(mcmc)) {
                stop(paste("mcmc is missing"))
            }
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(complete)) {
                stop(paste("Complete log-likelihood function is missing"))
            }
            if (missing(data)) {
                stop(paste("Data is missing"))
            }
            if (missing(sjwinit)) {
                sjwinit = 0
            }
            tpm <- proc.time()
            hold <- sjw(mcmc, z, complete, x = data, sjwinit,thrSJW,maxSJW)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]),3)
            #cat(paste("                             ", time, "  (converged at",  hold$iterations, "iterations)\n"))
	    myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,hold$status)
            timings[t] <- time
            t <- t + 1
        }


        if (method[l] == "DATA-BASED") {
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(data)) {
                stop(paste("Data is missing"))
            }
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- dataBased(x=data,K,z)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]),3)
            #cat(paste("                      ", time,"\n"))
	    myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,"OK")
            timings[t] <- time
            t <- t + 1
        }

        if (method[l] == "AIC") {
            if ((missing(constraint)==FALSE) && (length(constraint)>1)) {
		if (missing(mcmc)) {
			stop(paste("MCMC is missing"))
		}
                while (fic < length(constraint)) {
                  fic <- fic + 1
                  #if (fic == 1) {
                   # cat(paste(" (using constraint", fic, "of", length(constraint)),")", sep = "")
                  #}
                  #else {
                   # cat(paste("        (using constraint", fic, "of",length(constraint)), ")", sep = "")
                  #}
                  tpm <- proc.time()
                  permutations[[t]] <- aic(mcmc,constraint[fic])$permutations
                  time <- proc.time() - tpm
                  time <- round(as.numeric(time[3]),3) #as.numeric(time[3])
                  #cat(paste("   ", time, "\n"))
		  voutsas <- paste(method[l]," (constraint ", fic, " of ", length(constraint),")",sep="")
	          myPrettyPrint(gap0,gap1,gap2,gap3,voutsas,time,"OK")
                  timings[t] <- time
                  t <- t + 1
                }
            }
            else {
		    if (missing(mcmc)) {
		        stop(paste("MCMC is missing"))
		    }
		    if (missing(constraint)) {
		        constraint = 1 #default
		    }
		    tpm <- proc.time()
		    hold <- aic(mcmc,constraint)
		    permutations[[t]] <- hold$permutations
		    time <- proc.time() - tpm
		    time <- round(as.numeric(time[3]),3)
		    #cat(paste("                             ", time,"\n"))
       	            myPrettyPrint(gap0,gap1,gap2,gap3,method[l],time,"OK")
		    timings[t] <- time
		    t <- t + 1
                }
        }

       if (method[l] == "USER-PERM") {
	    #if(missing(userPerm)){cat('\n');stop(cat(paste("[ERROR]: user-defined permutations not defined"),'\n'))}
	    while (userIterator < userLength) {
		    tpm <- proc.time()
		    #if(userIterator==0){
			#cat(paste(" (input ", userIterator + 1, " of ", userLength,")", sep = ""))
		     #}else{
			#cat(paste("              (input ", userIterator + 1, " of ", userLength,")", sep = ""))
		    #}
		    permutations[[t]] <- userPerm[[userIterator+1]]
		    userPerm[[userIterator+1]] <- 0
		    time <- proc.time() - tpm
		    time <- round(as.numeric(time[3]),3)
		    if(userIterator>1){voutsas <- paste(method[l]," (input ", userIterator+1, " of ", userLength,")",sep="")}else{
			    voutsas <- paste(method[l])
		    }
		    if(userStatus[userIterator+1]==1){
	 	    	    myPrettyPrint(gap0,gap1,gap2,gap3,voutsas,"NA","OK")
		    }else{
			    myPrettyPrint(gap0,gap1,gap2,gap3,voutsas,"NA","FAIL")
		    }
		    #cat(paste("         NA","\n"))
		    timings[t] <- time
		    t <- t + 1
		    userIterator <- userIterator + 1
	    }
        }


    }
    cat(paste("    ......................................................................................\n"))
    cat(paste("\n"))

#    if ((dimname > 1) && (missing(z) == FALSE)){
    if ((missing(z) == FALSE)){
	if(missing(groundTruth)==TRUE){
		#by default the first method is chosen
		cat(paste("    Relabelling all methods according to method",nam[1], "..."))
		temp <- z
		for (i in 1:m) {
		    temp[i, ] <- order(permutations[[1]][i,])[z[i, ]]
		}
		zpivot <- apply(temp,2,function(y){uy <- unique(y);uy[which.max(tabulate(match(y, uy)))]}) # finds the mode
		myComparison <- compare.clust(zpivot, permutations, z, K)
		best.clusterings <- myComparison$clusters
		similarity.matrix <- myComparison$similarity[1:dimname,1:dimname]
		rownames(best.clusterings) <- nam
		#colnames(best.clusterings) <- nam
		permutations <- myComparison$permutations 
	}else{
		if( length(groundTruth) != dim(z)[2] ){
			stop(paste("groundTruth size not compatible"))
		}
		cat(paste("    Relabelling all methods according to ground truth ..."))
		temp <- z
		myComparison <- compare.clust(groundTruth, permutations, z, K)
		best.clusterings <- myComparison$clusters
		similarity.matrix <- myComparison$similarity
		rownames(best.clusterings) <- nam
		permutations <- myComparison$permutations 
	}
	results <- list(permutations, best.clusterings, timings, 
	similarity.matrix)
	names(results) <- c("permutations", "clusters", "timings", 
	"similarity")
    }#else{
#	results <- list(permutations,  timings)
#	names(results) <- c("permutations", "timings")
 #   }


    cat(paste(" done!\n"))
    cat(paste("    Retrieve the", dimname, "permutation arrays by typing:\n"))
    for (i in 1:dimname) {
        cat(paste("        [...]$permutations$\"", nam[i], "\"", sep = "", 
            "\n"))
    }
    cat(paste("    Retrieve the", dimname, "best clusterings: [...]$clusters\n"))
    cat(paste("    Retrieve the", dimname, "CPU times: [...]$timings\n"))
    moustakas <- 0;if(missing(groundTruth)==FALSE){moustakas <- 1}
    cat(paste("    Retrieve the", dimname + moustakas, "X", dimname+ moustakas, "similarity matrix: [...]$similarity\n"))
    runFinish <- proc.time() - runStart
    runFinish <- round(as.numeric(runFinish[3]),1) #as.numeric(time[3])
    cat(paste("    Label switching finished. Total time: ",runFinish," seconds.",sep=""),'\n')	
#    results <- list(permutations, best.clusterings, timings, 
#        similarity.matrix)
#    names(results) <- c("permutations", "clusters", "timings", 
#        "similarity")
    return(results)
}

compare.clust <- function (pivot.clust, perms, z, K){
    
    if (K < max(z)) {
        stop("K should be at least equal to max(z)")
    }
    if (dim(z)[2] != length(pivot.clust)) {
        if (K < max(z)) {
            stop("pivot.clust has not the same length with simulated z's")
        }
    }
    if(is.list(perms)==FALSE){
	l <- vector("list",length=1)
	l[[1]] <- perms
	names(l) <- "perm"
	perms <-l
    }
    outDim <- dim(summary(perms))[1]
    mySim <- numeric(outDim)
    ramone <- array(data = NA, dim =c(outDim,length(pivot.clust)))
    johnny <- perms
    for(joey in 1:outDim){
	n <- length(pivot.clust)
	m <- dim(z)[1]
	deedee <- array(data = 0,dim = c(m,n))
	k <- K
	st <- 1:k
	cost.matrix <- matrix(numeric(k * k), nrow = k, ncol = k)
	s <- 1:n
	for (iter in 1:m) {
		deedee[iter,] <- order(perms[[joey]][iter,])[z[iter, ]]
	}
	#ramone[joey,] <- apply(deedee,2,Mode)
	ramone[joey,] <- apply(deedee,2,function(y){uy <- unique(y);uy[which.max(tabulate(match(y, uy)))]}) # finds the mode
	alloc = ramone[joey,]
	for (i in 1:k) {
	    so <- pivot.clust[s[alloc == i]]
	    l <- length(alloc[alloc == i])
	    for (j in 1:k) cost.matrix[i, j] <- l - length(so[so == j])
	}
	matr <- lp.assign(cost.matrix)$solution
	perm <- numeric(k)
	for (i in 1:k) {
	    perm[i] <- st[matr[, i] > 0]
	}
	ramone[joey,] <- order(perm)[alloc]
	mySim[joey] <- length(which(ramone[joey,]==pivot.clust))
	for(i in 1:m){
		johnny[[joey]][i,] <- perms[[joey]][i,][perm]
	}
	#print(mySim[joey])
    }
    mySim<-mySim/m
    n <- length(pivot.clust)
    similarity.matrix <- array(data = 0, dim = c(outDim+1, outDim+1))
    for (i in 1:(outDim)) {
        for (j in 1:i) {
            similarity.matrix[i, j] <- length(which(ramone[i,]==ramone[j,]))
        }
    }
    i <- outDim+1
    for (j in 1:(outDim)) {
        similarity.matrix[i, j] <- length(which(pivot.clust == ramone[j,]))
    }
    similarity.matrix = (similarity.matrix + t(similarity.matrix))/n
    diag(similarity.matrix) <- rep(1, outDim+1)
    rownames(similarity.matrix) <- c(names(perms), deparse(substitute(pivot.clust)))
    colnames(similarity.matrix) <- c(names(perms), deparse(substitute(pivot.clust)))
    results <- list(similarity.matrix, ramone,johnny)
    names(results) <- c("similarity", "clusters","permutations")
    return(results)
}

