`fblr.weight` <-
function(y,bin,niter,thin=5,nburn=10000, kmax = 10, geo=1, 
		delta1=0.001,delta2= 0.1, predict=FALSE, group=NULL,weight = NULL, 
		file ="fblr_mcmc.txt"){

	require(MASS) # library MASS for mvrnorm

	# Checks
	if(!is.null(group)){
   		if(!all(names(table(unlist(group))) == 1:ncol(bin)))
    			stop("Not every binary contained in group")
    	}

	if(!is.null(weight)){
  		if(length(weight) != ncol(bin))
  			stop("length of weight and number of binaries differs")
  	}

	# Auxiliary functions

	In <- function(n) diag(rep(1,n))

	int.matrix <- function(group){
    		nbin <- length(unique(unlist(group)))
    		IM <- matrix(0,ncol=nbin,nrow=nbin)
    		for(i in 1:length(group)){
        		ingroup <- (1:nbin %in% group[[i]])*1
        		IM.group <- ingroup%*%t(ingroup)
        		IM <- IM + IM.group
        	}
        	IM <- (IM > 0)*1
        	diag(IM) <- 0
		IM
	}

	weight.matrix <- function(wst,IM){
    		WM <- outer(wst,wst,FUN=function(x,y) (x+y)/2 )
    		WIM <- WM * IM   # only existing interactions
    	WIM
	}

	# Generate z values from truncated normal with the Inverse cdf method
	trunc.norm <- function(mu){
    		mu + qnorm( (pnorm(-mu) +runif(length(mu))*(1-pnorm(-mu))))
    	}
	generate.z <- function(y,eta){
    		ifelse(y == 1, trunc.norm(eta), -trunc.norm(-eta))
	}

	# the non-Bayesfactor part of the acceptance probability
  	compute.R <- function(model, mi, geo){
        	R <- 1
        	sizes <- model[mi$begin]
        	if (model[1] == 0 | mi$model[1] == 0) { # move to or from Null Model
            		R <- switch(2 + 1 * (model[1] == 0) - 1 * (mi$model[1] == 
                		0), 4, 1, 1/4)
        	}
        
        	else if (model[1] == mi$model[1]){     # change
            		if (sum(sizes) < sum(mi$model[mi$begin])) 
                		R <- geo                       # cbirth
            		else if (sum(sizes) > sum(mi$model[mi$begin])) 
                		R <- 1/geo                     # cdeath
        	}
        	else if (model[1] < mi$model[1]) 
            		R <- 1/(sum(sizes == 1) + 1)       # birth
        	else if (model[1] > mi$model[1]) 
            		R <- sum(sizes==1)                 # death
        	R
    	}

	#Compute acceptance probability
	acceptance5 <- function(model.old, move.info ,marg.loglike,marg.loglike.new,geo){
    		min(c(1, compute.R(model.old, move.info,geo=geo)*  move.info$w *
    			exp(marg.loglike.new -marg.loglike)))# Bayes Factor
    	}


	compare <- function(x,y) (all(x == y))*1

	eval.logic.single <- function(logic,bin){
    		size <- logic[1]
    		if(size == 0) return(NULL)
    		lvars <- logic[2:(size+1)]
    		norm <- (lvars > 0)*1       # norm means that variables are not negated
    		if(size == 1)  
			return( (bin[,abs(lvars)] == norm)*1)
    		else 
			apply(bin[,abs(lvars)],1,compare,y = norm)
    	}

	# update model matrix
	change.matrix <- function(mi, B) {
        	if(mi$model[1] == 0) 
            		return(matrix(rep(1, n), ncol = 1))
        	if(mi$move.type == 3) 
            		B[, mi$pick + 1] <- eval.logic.single(mi$model[mi$begin[mi$pick]:(mi$begin[mi$pick] + 2)], bin = bin)
        	if(mi$move.type == 2) 
            		B <- cbind(B, if (mi$k == 0) 
                		eval.logic.single(mi$model[2:4], bin = bin)
            			else eval.logic.single(mi$model[(mi$begin[mi$k] + 3):
					(mi$begin[mi$k] + 5)], bin = bin))
        	if(mi$move.type == 1) {
            		B[, mi$pick + 1] <- B[, mi$k + 1]
            		B <- B[, 1:mi$k]
        	}
        	B
    	}

    	# generate a candidate model
	gen.cand <- function(model){ 
        	k <- model[1]
        	if(k == 0){
                	lvar <- sample(vars,1)
                	model[1:3] <- c(1, 1, sample(c(-1, 1), 1) * lvar)
                	return(list(model=model,k=0, begin=NULL,  move.type=2, pick=NULL, 
				w=wst[lvar]))    
                }
        	begin <- seq(from = 2, length.out = k, by = 3)
        	sizes <- model[begin]
        	death.poss <- sum(sizes == 1)
            	u <- runif(1)
            	if (u < (death.poss > 0) * 0.25) 
                	move.type <- 1
            	else if (u < (death.poss > 0) * 0.25 + (k < kmax) * 0.25) 
                	move.type <- 2
            	else move.type <- 3
      		switch(move.type, death(model,k,begin,sizes,death.poss), birth(model,k,begin), 
                	change(model, k, begin, sizes))  
    	}

 	# Functions for the different moves
 	death <- function(model,k, begin, sizes, death.poss){
        	dead <- ifelse(death.poss == 1, which(sizes == 1), sample(which(sizes == 1), 1))
        	w <- 1/wst[abs(model[begin[dead]+1])]
        	model[begin[dead]:(begin[dead] + 2)] <- model[begin[k]:(begin[k] + 2)]
        	model[begin[k]:(begin[k] + 2)] <- rep(0,3)
        	model[1] <- model[1] - 1
        	list(model=model, k=k,begin=begin, move.type=1, pick= dead,w=w)
    	}

	birth <- function(model,k,begin){
        	born <- sample(vars,1)
        	model[(begin[k] + 3):(begin[k] + 4)] <- c(1, sample(c(-1, 1), 1) * born)
        	model[1] <- model[1] + 1
        	list(model=model,k=k,begin=begin, move.type=2, pick=NULL,w = wst[born])
    	}

	change <- function(model,k,begin,sizes){
        	pick <- sample(1:k, 1)
        	size <- sizes[pick]
        	u <- runif(1)
        	if (u < 0.5* (size > 1)) cmove <- 1
        	else if (u < 0.5 * (size != 2)) cmove <- 2
        	else cmove <- 3
        	logic <- model[begin[pick]:(begin[pick] + 2)] 
        	change.info <- switch(cmove, cdeath(logic), cbirth(logic), cchange(logic))  
        	model[begin[pick]:(begin[pick] + 2)] <- change.info$logic
        	list(model=model, k=k,begin=begin,move.type=3,pick=pick,w=change.info$w)
    	}



	cdeath <- function(logic) {
        	size <- logic[1]
        	wij <- WIM[abs(logic[2]),abs(logic[3])] # weight of the two-factor-interaction
        	pickvar <- sample(2:3) # pickvar[1] is thrown out, pickvar[2] stays in the model
        	ch <- abs(logic[pickvar[2]])
        	logic[c(1, pickvar[1], size + 1)] <- c(size - 1, logic[size + 1], 0)
        	w <- wst[ch]/(wij * ints[ch])
        	list(logic=c(logic[1], sort.abs(logic[2:size]), rep(0, 2 - size + 1)), w=w)
    	}

	cbirth <- function(logic){
     		size <- logic[1]
     		lvar <- abs(logic[2])
     		if(sum(IM[lvar,])==0) return(cchange(logic))
     		newvar <- ifelse(sum(IM[lvar,]> 0)== 1, vars[IM[lvar,]>0], 
			sample(vars[IM[lvar,]>0],1))
     		# No sampling, if there is only one possibility
     		w <- ints[lvar] * WIM[lvar,newvar] / wst[lvar]
     		list(logic=c(size+1,sort.abs(c(logic[2],sample(c(-1,1),1)*newvar))), w=w)
     	}


	cchange <- function(logic) {
    		size <- logic[1]
    		if(size == 1){
        		newvar <- sample(vars,1)
        		w <- wst[newvar]/wst[abs(logic[2])]
    			return(list(logic=c(1,sample(c(-1,1),1)*newvar,0),w=w))
    		} 
		else{
        		pickvar <- logic[sample(2:3)]  # pickvar[1] is thrown out, pickvar[2] stays in the model
        		newvar <- ifelse(sum(IM[abs(pickvar[2]),])==1, 
				vars[IM[abs(pickvar[2]),]>0], sample(vars[IM[abs(pickvar[2]),]>0],1))
        		w <- WIM[abs(pickvar[2]),newvar] / WIM[abs(pickvar[1]),abs(pickvar[2])]
    			return(list(logic=c(2,sort.abs(c(pickvar[2], sample(c(-1,1),1)*newvar))),w=w))
          	}
	}

	sort.abs <- function(vek) vek[order(abs(vek))]


	write.table(file=file, x=NULL, col.names=FALSE) # Delete old file

    	# Initialize
    	nbin <- dim(bin)[2]
    	vars <- 1:dim(bin)[2]
    	n <- length(y)
    	if(is.null(group)) group <- list(1:nbin)
    	if(is.null(weight)) weight <- rep(1,nbin)
    	IM <- int.matrix(group)
    	wst <- weight/sum(weight)*length(weight) # scale weights
    	WIM <- weight.matrix(wst,IM)
    	ints <- apply(IM,1,function(x) sum(x > 0))
    	ints <- ints/mean(ints)
    	if(predict) pred <- rep(0,n)
    	# Initialize
    	model <- rep(0,1+kmax*(3))
    	k <- 0
    	B <- as.matrix(rep(1,n))
   	b <- 0
    	accept <- 0
    	# MCMC
    	for(i in 1:(niter+nburn)){
    		# Update z
    		eta <- B%*%b
    		z <- generate.z(y,eta)
    		# Update tau
    		tau <- rgamma(1,shape = delta1 + 0.5*(k+1), rate = delta2 + 0.5* t(b)%*%b)
    		Vstar <- solve(t(B)%*%B + tau*In(k+1))
    		bstar <-  t(z)%*%z - t(z)%*%B%*%Vstar%*%t(B)%*%z
    		marg.loglike <- 0.5*(log(abs(det(Vstar))) - bstar + (k+1)*log(tau))
   	 
    		# Generate new model
    		move.info <- gen.cand(model)
    		B.new <- change.matrix(move.info,B=B)
    		k.new <- dim(B.new)[2]-1
    		Vstar.new <-  solve(t(B.new)%*%B.new + tau*In(k.new+1))
    		bstar.new <-  t(z)%*%z - t(z)%*%B.new%*%Vstar.new%*%t(B.new)%*%z
    		marg.loglike.new <- 0.5*(log(abs(det(Vstar.new))) - bstar.new + (k.new+1)*log(tau))
    		if(runif(1) < acceptance5(model, move.info, marg.loglike,marg.loglike.new,geo=geo)){
    			model <- move.info$model
    			B <- B.new
    			k <- k.new
    			Vstar <- Vstar.new
    			marg.loglike <- marg.loglike.new
     			accept <- accept + 1
    		}
     		# Update beta
    		m <-  Vstar%*%t(B)%*%z
    		b <- mvrnorm(1,mu=m,Sigma=Vstar)
		if(i%%thin == 0 & i > nburn){
    			write.table(matrix(c(model,round(marg.loglike,2),round(tau,6),
				round(b,4),rep(0,kmax+1-length(b))),nrow=1),file=file,
    				append=TRUE,row.names=FALSE,col.names=FALSE)
    			if(predict) pred <- pred + pnorm(B%*%b)
    		}
	}
	if(predict) 
		return(list(accept=accept/(niter+nburn), pred=as.vector(pred)/(niter/thin)))
	else
		return(accept/(niter+nburn))
}

