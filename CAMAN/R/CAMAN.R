mixalg = function(obs, weights=NULL, family="gaussian", data=NULL, 
                  pop.at.risk=NULL, var.lnOR=NULL, limit=0.01, 
                  acc=10^(-7), numiter=5000, startk=50){        
    # Performs CAMAN (computer-assisted analysis of mixtures)
    if (family == "gaussian") dens_i = 0
    else if (family == "poisson") dens_i = 1
    else if (family == "binomial") dens_i = 2
    else return("Please enter a valid family of distribution (normal, poisson, binomial)")    
    
    #check data
    if (is.null(data)) data <- data.frame() # no data was given 
    if (!((obs %in% colnames(data))||is.numeric(obs))) stop("obs must be a colname of 'data' or a numeric vector")
    if (!((weights %in% colnames(data))||is.numeric(weights) || is.null(weights))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
    if (!((var.lnOR %in% colnames(data))||is.numeric(var.lnOR) || is.null(var.lnOR))) stop("variances must be a colname of 'data', a numeric vector or 'NULL'")

	if (is.null(var.lnOR) ) is_metaAnalysis <- 0 #if variances are not specified, we apply don't perform a meta analysis and estimate the variances
	else is_metaAnalysis <- 1 #variances are given --> meta analysis --> no estimation

    
    n= max(nrow(data), length(obs) )
    datmat <- matrix(1,ncol=4, nrow=n)
    
    if (is.numeric(obs)&&length(obs>1)) datmat[,1] <- obs
    else datmat[,1] <- data[,obs]
    
    #build matrix 'datmat' by reading out the command
    tmpdat <- list(obs, weights, pop.at.risk, var.lnOR)
    for (i in 1:4){
        if (is.character(tmpdat[[i]]) ) datmat[,i] <- data[,tmpdat[[i]]] #colname was given
        else if (is.null(tmpdat[[i]]) ) datmat[,i] <- rep(1,n)   #NULL was given
        else if (is.numeric(tmpdat[[i]]) ) datmat[,i] <- tmpdat[[i]] #a numeric vector was given
        else stop("Data initialization failed...")
    }
    rm(tmpdat)
  #estimate variances  
	if ((sum(rep(1,n) == datmat[,4])== n) && family=="gaussian" && is.null(var.lnOR)) datmat[,4] <- rep(var(datmat[,1]),n)
    
	
    res1 <- .C("caman_C", as.vector(as.double(datmat[,1])), as.vector(as.double(datmat[,2])), as.vector(as.double(datmat[,3])), 
        as.vector(as.double(datmat[,4])), as.integer(n), as.integer(startk), as.integer(dens_i), 
        as.integer(999), as.double(999), rep(as.double(999), 150), rep(as.double(999), 150), 
        as.double(limit), as.double(acc), as.integer(numiter), as.double(c(-999)), 
               rep(as.double(-999), (2*startk + 2)),rep(as.double(-999), 2) ,
               as.integer(is_metaAnalysis)  ,PACKAGE = "CAMAN")
    
    numObs <- sum(datmat[,2])
    k <- res1[[8]]
    p=res1[[10]][1:k]    
    bic <- -2 * res1[[9]] + (2*k - 1) * log(numObs)
    VEM_tmp <- res1[[16]]
    EM_tmp <- res1[[17]]
    finalacc <- c(VEM_tmp[2], EM_tmp[1]) # VEM, EM 
    vem_res <- matrix(VEM_tmp[3: (2*VEM_tmp[1] + 2)], ncol=2)
    totalsteps <- c(res1[[14]][1], EM_tmp[2]) #VEM, EM
    res <- new("CAMAN.object",dat=datmat, family=family, LL=res1[[9]], 
               num.k=k, p=p, t=res1[[11]][1:k], num.obs = numObs, steps=totalsteps, 
               otherParams = c(limit, numiter, acc, startk), BIC = bic, VEM_result = 
                 vem_res, finalacc= finalacc, is_metaAnalysis=is_metaAnalysis)
    if (dens_i == 0) res@component.var=res1[[15]]
    
	#compute posterior prob
	probs <- mix.densDistr(res)
    res@prob <- probs
    res@classification <- apply(probs, 1, which.max)
    
    
    if (res@steps[1] >= res@otherParams[2]) {
      warning("Warning: Solution does not satisfy convergence criterion:\n The last VEM-iteration had a accuracy of ",
              res@finalacc[1],". You asked for (acc=)", acc,
              "\n Please increase numiter or decrease acc", sep="")}
if (res@steps[2] >= res@otherParams[2]) {
  warning("Warning: Solution does not satisfy convergence criterion:\n The last EM-iteration had a accuracy of ",
          res@finalacc[2],". You asked for (acc=)", acc,
          "\n Please increase numiter or decrease acc", sep="")}
    return(res)
}


anova.CAMAN.object <- function(object, object1, nboot=2500, limit=0.01, 
                               acc=10^(-7), numiter=5000, 
                               giveBootstrapData=FALSE, giveLikelihood=FALSE, ...){ 
mix0 <- object
mix1 <- object1
#compute LL-ratio:
# simulate data from mix0 and compare LL of mix0 and mix1 based on this data   
	cl <- match.call()	
    if (nboot<1) return("Please enter a valid number for nboot")
    
    if (mix0@family == "gaussian") dens_i0 = 0
        else if (mix0@family == "poisson") dens_i0 = 1
        else if (mix0@family == "binomial") dens_i0 = 2

    if (mix1@family == "gaussian") dens_i1 = 0
        else if (mix1@family == "poisson") dens_i1 = 1
        else if (mix1@family == "binomial") dens_i1 = 2
        
    datUnpack <- rep(mix0@dat[,1], mix0@dat[,2])
    len <- length(datUnpack)
    obs <- sapply(1:nboot, function(x){simDat(mix0)}) #returns a matrix with bootstrapped observations of m1
	
    #do we need to unpack the data...
    if (sum(mix0@dat[,2]) == nrow(mix0@dat) ){ #no weights --> datUnpacked == mix@dat[,1]
		dataUnpacked=FALSE
    }
    else {  #else: datUnpacked != mix@dat[,1] --> there were weights
		dataUnpacked=TRUE
#if the got unpacked before, we now pack them together again
		tmp_obs <- apply(obs,2, function(x){as.numeric(names(table(x)))})
		tmp_weights <- apply(obs,2, function(x){as.numeric(table(x))})
		SimulatedObs <- tmp_obs
		SimulatedWeights <- tmp_weights
	}	
	
    LL0 <- NULL
    LL1 <- NULL
    for (i in 1:nboot){    
    #perform EM for each bootstrap sample
    #perform EM for mix0
		if (dataUnpacked){ #if the data was packed, we 
#need to handle it in another way because there might be less rows in our datamatrix... 
			obs_i <- SimulatedObs[[i]]
			obs_weights <- SimulatedWeights [[i]]
			tmplen <- length(obs_weights)
			col3 <- rep(1, tmplen) #packing--> use ones as parameter
			col4 <- rep(1, tmplen) #packing--> use ones as parameter
		}
		else { #no packing--> use original parameters for var.lnOR and pop.at.risk
			obs_i <- obs[,i]
			obs_weights <- mix0@dat[,2] #=rep(1, len)
			col3 <- mix0@dat[,3]
			col4 <- mix0@dat[,4]
			tmplen <- len			
		}
        res0 <- .C("mixalg_sub", as.double(obs_i), as.double(obs_weights), as.double(col3), 
            as.double(col4), as.integer(tmplen), as.integer(mix0@num.k), as.integer(dens_i0), 
            as.integer(mix0@num.k), as.double(999), as.double(mix0@p), as.double(mix0@t), as.double(limit), as.double(acc), 
            as.integer(numiter), as.double(c(-999)), as.integer(1),  as.integer(mix0@is_metaAnalysis) , PACKAGE = "CAMAN")   
        res1 <- .C("mixalg_sub", as.double(obs_i), as.double(obs_weights), as.double(col3), 
            as.double(col4), as.integer(tmplen), as.integer(mix1@num.k), as.integer(dens_i1), 
            as.integer(mix1@num.k), as.double(999), as.double(mix1@p), as.double(mix1@t), as.double(limit), as.double(acc), 
            as.integer(numiter), as.double(c(-999)), as.integer(1) , as.integer(mix1@is_metaAnalysis), PACKAGE = "CAMAN")          

		LL0[i] <- res0[[9]]
        LL1[i] <- res1[[9]]
    }
res <- list() 
LL_ratios <- sort(-2*(LL0 - LL1))  #90, 95, 97.5, 99 - quartils
LL_ratio_quartils <- LL_ratios[floor(c(.9, .95, .975, .99)*nboot)]
names(LL_ratio_quartils) <- c(.9, .95, .975, .99)
res$overview <- data.frame(c(as.character(cl$object), as.character(cl$object1)), 
                           c(mix0@num.k, mix1@num.k), c(mix0@BIC, mix1@BIC), 
                           c(mix0@LL, mix1@LL), c(NA, -2*(mix0@LL - mix1@LL)))
names(res$overview) = c("mixture model","k","BIC","LL", "LL-ratio")
res$"LL ratios in bootstrap-data" <- LL_ratio_quartils
res$"simulated p-value" <- sum(LL_ratios > (-2*(mix0@LL - mix1@LL)) )/nboot

if (giveBootstrapData) {
	if (dataUnpacked) res$BootStrapData <- SimulatedObs
	else res$BootStrapData <- obs
}
if (giveLikelihood) res$LL <- rbind(LL0,LL1) 

return(res)
}

mixalg.paraBoot <- function(mix.estim, nboot=50, limit=0.01, acc=10^(-7), numiter=5000, startk=50, giveBootstrapData=FALSE){
    # Performs a parametric bootstrap on data.
    # Returns the standard deviation of patameters t and p of the given mixture model! 
    if (nboot<1) return("Please enter a valid number for nboot")
    
    if (mix.estim@family == "gaussian") dens_i = 0
        else if (mix.estim@family == "poisson") dens_i = 1
        else if (mix.estim@family == "binomial") dens_i = 2

    datUnpack <- rep(mix.estim@dat[,1], mix.estim@dat[,2])
    len <- length(datUnpack)
    obs <- sapply(1:nboot, function(x){simDat(mix.estim)}) #returns a matrix with bootstrapped observations
    
    #other parameters for the EM-algorithm
    if (sum(mix.estim@dat[,2]) == nrow(mix.estim@dat) ){ #no weights --> datUnpacked == mix@dat[,1] ==> use original parameters for var.lnOR and pop.at.risk
        c2 <- mix.estim@dat[,2]
        c3 <- mix.estim@dat[,3]
        c4 <- mix.estim@dat[,4] 
    }
    else {  #else: there were weights --> datUnpacked != mix@dat[,1] ==> use ones as add. parameters!
        c2 <- rep(1, len)
        c3 <- c2
        c4 <- c2           
    }

    p_mat <- matrix(0, ncol=mix.estim@num.k, nrow=nboot)
    t_mat <- matrix(0, ncol=mix.estim@num.k, nrow=nboot)

cat("\nProgress:\n0.........50.......100% of Bootstraps done.\n")
    for (i in 1:nboot){    
    #perform EM for each bootstrap sample
        res1 <- .C("mixalg_sub", as.double(obs[,i]), as.double(c2), as.double(c3), 
        as.double(c4), as.integer(len), as.integer(mix.estim@num.k), as.integer(dens_i), 
        as.integer(mix.estim@num.k), as.double(999), as.double(mix.estim@p), as.double(mix.estim@t), as.double(limit), as.double(acc), 
        as.integer(numiter), as.double(c(-999)), as.integer(1) , as.integer(mix.estim@is_metaAnalysis), PACKAGE = "CAMAN")   
    if (i %in% seq(0,nboot,max(floor(nboot/20),1))  ) cat ("|")

    p_mat[i, ] <- res1[[10]][1:mix.estim@num.k]
    t_mat[i, ] <- res1[[11]][1:mix.estim@num.k]
    }
cat ("\n")
return(list(p_mat, t_mat, obs=obs))
    sd.p <- apply(p_mat, 2, sd)
    sd.t <- apply(t_mat, 2, sd)
    if (giveBootstrapData) return(list(sd.p = sd.p, sd.t = sd.t, giveBootstrapData = obs))
    return(list(sd.p = sd.p, sd.t = sd.t))
}


mixalg.boot <- function(mix, nboot=500, limit=0.01, acc=10^(-5), numiter=5000, startk=50, returnBootstrapRep= FALSE){
    # Performs a nonparametric bootstrap on data.
    # Returns the computed optimal number of component for each bootstrap replication 
    if (nboot<1) return("Please enter a valid number for nboot")
    
    if (mix@family == "gaussian") dens_i = 0
    else if (mix@family == "poisson") dens_i = 1
    else if (mix@family == "binomial") dens_i = 2

    datUnpack <- rep(mix@dat[,1], mix@dat[,2]) # If there are weights --> unpack data and do not use any weights/ -- freq = 1
    #generate bootstrap samples:

    x.sample <- matrix( rep(datUnpack,nboot), ncol=nboot, byrow=FALSE) #initialization
    xlen <- nrow(x.sample)
    x.sample <- apply(x.sample, 2, function(yy){sample(yy, xlen, replace=TRUE)}) #sample (bootstrap)
	
	#pack data
    tmp_obs <- apply(x.sample,2, function(x){as.numeric(names(table(x)))})
    tmp_weights <- apply(x.sample,2, function(x){as.numeric(table(x))})
    vec_n <- sapply(tmp_obs, length)
    bootSamples <- unlist(tmp_obs)
    bootWeights <- unlist(tmp_weights)
    
    tmpn <- nrow(mix@dat)
    bootVar <- rep(1, length(bootSamples))
    if (sum(mix@dat[,4] == rep(1,tmpn)) != tmpn)
    {  #variances were given --> extract them
        tmpvar <- mix@dat[,4]
        names(tmpvar) <- mix@dat[,1] 
        tmpBootVar <- sapply(tmp_obs, function(x){tmpvar[as.character(x)]})
        bootVar <-  as.numeric(unlist(tmpBootVar))
    }
    bootPopAtRisk <- rep(1, length(bootSamples))
    if (sum(mix@dat[,3] == rep(1,tmpn)) != tmpn)
    {  #popAtRisk were given --> extract them
        tmppop <- mix@dat[,3]
        names(tmppop) <- mix@dat[,1] 
        tmpBootPop <- sapply(tmp_obs, function(x){tmppop[as.character(x)]})
        bootPopAtRisk <- as.numeric(unlist(tmpBootPop))
    }
    
    #permutation step; write observations in a vector (columnwise)    
    res1 <- .C("caman_boot", as.double(bootSamples), as.double(bootWeights), as.double(bootPopAtRisk), 
        as.double(bootVar), as.vector(as.integer(vec_n)), as.integer(startk), as.integer(dens_i), 
        as.integer(999), rep(as.double(999.99), nboot), rep(as.double(999), 150), rep(as.double(999), 150), 
        as.double(limit), as.double(acc), as.integer(numiter), as.double(c(-999)), as.integer(nboot), rep(as.integer(999),nboot), rep(as.double(999.99), nboot), as.integer(mix@is_metaAnalysis),PACKAGE = "CAMAN")
    
    if (returnBootstrapRep) res <- list(dat.bootstrap=x.sample, LL=res1[[9]], numk.boot=res1[[17]], LL_k1 = res1[[18]]) #unpacked bootstrap data is returned! 
    else res <- list(LL=res1[[9]], numk.boot=res1[[17]], LL_k1 = res1[[18]])
    return(res)
}




mixalg.EM <- function(mix = NULL, p, t, obs=NULL, weights=NULL, family="gaussian", 
                      data=NULL, pop.at.risk=NULL, var.lnOR=NULL,  limit=0.01, 
                      acc=10^(-7), numiter=5000){
    # computes the seconde (EM-) part of the CAMAN algorithm: 
    # use manualy defined values for p, t &  k and 
    #   --> refiend soultion with EM algorithm 
    # return updated estimates for p, t, acc & number of iterations (numiter)
    if (length(p) == length(t) ) num.k = length(t)
    else stop("Please enter valid data for p and t")
    if (!is.null(mix)){
        family = mix@family
        datmat = mix@dat   
		is_metaAnalysis = mix@is_metaAnalysis
    }
    else{     
	    #check data
	    if (is.null(data)) data <- data.frame() # no data was given 
	    if (!((obs %in% colnames(data))||is.numeric(obs))) stop("obs must be a colname of 'data' or a numeric vector")
	    if (!((weights %in% colnames(data))||is.numeric(weights) || is.null(weights))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
	    if (!((var.lnOR %in% colnames(data))||is.numeric(var.lnOR) || is.null(var.lnOR))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
 
 	    if (is.null(var.lnOR) ) is_metaAnalysis <- 0 #if variances are not specified, we apply don't perform a meta analysis and estimate the variances
    	else is_metaAnalysis <- 1 #variances are given --> meta analysis --> no estimation

	    
	    n= max(nrow(data), length(obs) )
	    datmat <- matrix(1,ncol=4, nrow=n)
	    
	    if (is.numeric(obs)&&length(obs>1)) datmat[,1] <- obs
	    else datmat[,1] <- data[,obs]
		
		
	    
	    #build matrix 'datmat' by reading out the command
	    tmpdat <- list(obs, weights, pop.at.risk, var.lnOR)
	    for (i in 1:4){
	        if (is.character(tmpdat[[i]]) ) datmat[,i] <- data[,tmpdat[[i]]] #colname was given
	        else if (is.null(tmpdat[[i]]) ) datmat[,i] <- rep(1,n)   #NULL was given
	        else if (is.numeric(tmpdat[[i]]) ) datmat[,i] <- tmpdat[[i]] #a numeric vector was given
	        else stop("Data initialization failed...")
	    }
		#estimate variances  
		if ((sum(rep(1,n) == datmat[,4])== n) && family=="gaussian" && is.null(var.lnOR)) datmat[,4] <- rep(var(datmat[,1]),n)
		
		rm(tmpdat)
    }
	if (family == "gaussian") dens_i = 0
	else if (family == "poisson") dens_i = 1
	else if (family == "binomial") dens_i = 2    
	else stop("Please enter a valid density distribution (gaussian, poisson, binomial)")
	
	
    res1 <- .C("mixalg_sub", as.double(datmat[,1]), as.double(datmat[,2]), as.double(datmat[,3]), 
    as.double(datmat[,4]), as.integer(nrow(datmat)), as.integer(num.k), as.integer(dens_i), 
    as.integer(num.k), as.double(999), as.double(p), as.double(t), as.double(limit), as.double(acc), 
    as.integer(numiter), as.double(c(-999)), as.integer(1), as.integer(is_metaAnalysis) ,PACKAGE = "CAMAN")
    
    numObs <- sum(datmat[,2])
    bic <- -2 * res1[[9]] + (2*num.k - 1) * log(numObs)
    totalsteps <- c(NA, res1[[14]]) #VEM, EM
    finalacc <- c(NA, res1[[13]])
    
    res <- new("CAMAN.object",dat=datmat, family=family, LL=res1[[9]], num.k=num.k, p=res1[[10]], t=res1[[11]], 
			num.obs = numObs, steps=totalsteps, otherParams = c(limit, numiter, acc, startk=num.k), BIC = bic, 
			VEM_result = matrix(), finalacc= finalacc, is_metaAnalysis = is_metaAnalysis)
    if (dens_i == 0) {
		if  (num.k>1) res@component.var=res1[[ 15 ]]
		else res@component.var = var(res@dat[,1])
	}
	
    #compute posterior probabilities
    probs <- mix.densDistr(res)
    res@prob <- probs
    res@classification <- apply(probs, 1, which.max)
    
return(res)
}


mixalg.VEM <- function(mix = NULL, obs=NULL, weights=NULL, data=NULL, pop.at.risk=NULL, var.lnOR=NULL, family="gaussian", limit=0.01, acc=10^(-7), numiter=5000, startk=50){
    # computes the first part of the CAMAN algorithm:
    #    1. construct grid of potential subpopulation means 
    #    2. calculation of mxing kernel density
    #    3. VEM algorithm 
    # returns estimates for parameters t & weights p
	
    if (!is.null(mix)){
		family = mix@family
		datmat = mix@dat   
		is_metaAnalysis = mix@is_metaAnalysis
	}
	else{     
		#check data
		if (is.null(data)) data <- data.frame() # no data was given 
		if (!((obs %in% colnames(data))||is.numeric(obs))) stop("obs must be a colname of 'data' or a numeric vector")
		if (!((weights %in% colnames(data))||is.numeric(weights) || is.null(weights))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
		if (!((var.lnOR %in% colnames(data))||is.numeric(var.lnOR) || is.null(var.lnOR))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
		
		if (is.null(var.lnOR) ) {
			is_metaAnalysis <- 0 #if variances are not specified, we apply don't perform a meta analysis and estimate the variances
		}
		else {
			is_metaAnalysis <- 1 #variances are given --> meta analysis --> no estimation
		}
		
		n = max(nrow(data), length(obs) )
		datmat <- matrix(1,ncol=4, nrow=n)
			
		if (is.numeric(obs)&&length(obs>1)) datmat[,1] <- obs
		else datmat[,1] <- data[,obs]
		
		#build matrix 'datmat' by reading out the command
		tmpdat <- list(obs, weights, pop.at.risk, var.lnOR)
		for (i in 1:4){
			if (is.character(tmpdat[[i]]) ) datmat[,i] <- data[,tmpdat[[i]]] #colname was given
			else if (is.null(tmpdat[[i]]) ) datmat[,i] <- rep(1,n)   #NULL was given
			else if (is.numeric(tmpdat[[i]]) ) datmat[,i] <- tmpdat[[i]] #a numeric vector was given
			else stop("Data initialization failed...")
		}
		#estimate variances  
		if ((sum(rep(1,n) == datmat[,4])== n) && family=="gaussian" && is.null(var.lnOR)) datmat[,4] <- rep(var(datmat[,1]),n)
		
		rm(tmpdat)
	}
    if (family == "gaussian") dens_i = 0
    else if (family == "poisson") dens_i = 1
    else if (family == "binomial") dens_i = 2
    else stop("Please enter a valid density distribution (normal, poisson, binomial)")

    num.k <- startk #min(sum(datmat[,2]),startk)
    p <- rep(999, num.k)
    t <- rep(999, num.k)
	#perform the EM algorithm for the given (simulated) data
    res1 <- .C("mixalg_sub", as.double(datmat[,1]), as.double(datmat[,2]), as.double(datmat[,3]), 
    as.double(datmat[,4]), as.integer(nrow(datmat)), as.integer(num.k), as.integer(dens_i), 
    as.integer(num.k), as.double(999), as.double(p), as.double(t), as.double(limit), as.double(acc), 
    as.integer(numiter), as.double(c(-999)), as.integer(0) , as.integer(is_metaAnalysis), 
	as.double(rep(-999.9,num.k)),PACKAGE = "CAMAN")

    LL <- res1[[9]]
	numObs <- sum(datmat[,2])
    bic <- -2 * res1[[9]] + (2*num.k - 1) * log(numObs)
    grid <- data.frame(p=res1[[10]], t=res1[[11]])
	grid <- grid[grid$p>0,]
	rownames(grid) <- as.character(1:nrow(grid))
	finalacc <- res1[[13]][1]
	totalsteps <- res1[[14]][1]
    
	totalgrid <- data.frame(p=res1[[10]], t=res1[[11]], gradient=res1[[18]])

	res <- new("CAMAN.VEM.object",dat=datmat, family=family, LL=LL, 
			num.k=num.k, num.obs = numObs, steps=totalsteps, otherParams = c(limit, numiter, acc, startk), 
			BIC = bic, finalacc= finalacc, startk=startk, grid=grid, totalgrid=totalgrid)
	
	
	
	return(res)    
}



mix.densDistr <- function(mix){
   # computes the probability for each observation (1..n - row of mix@dat) belonging to each component (1..k)
   # returns a matrix of dimension n x k
   dat <- mix@dat[,1]
   res <- matrix(ncol=mix@num.k, nrow=length(dat))
   p <- mix@p
  if (mix@family == "gaussian") {
       mu <- mix@t
       mix.sd <- sqrt(mix@component.var)
       for (i in 1:mix@num.k) res[,i] <- sapply(dat,
			function(x){p[i]*dnorm(x, mu[i], mix.sd ) / sum(p*dnorm(x, mu,
			mix.sd ))})
       }
  if (mix@family == "binomial") {
       prob <- mix@t
       popAtRisk <- mix@dat[,3]
       for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1,
			function(x){p[i]*dbinom(x[1], x[2], prob[i]) / sum(p*dbinom(x[1],
			x[2], prob))})
       }
  if (mix@family == "poisson") {
       lambda <-  mix@t
       popAtRisk <- mix@dat[,3]
       for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1, 
					function(x){p[i]*dpois(x[1], x[2]* lambda[i]) / 
								sum(p*dpois(x[1], x[2] * lambda))})
       }
   return(res)
}


getFDR <- function(dat, threshold=.7, idxNotDiff=1 ){
   #computes False Discovery Rate, etc. 
   tau0 <- dat@prob[,idxNotDiff] #p(not differentail genes)
   tau1 <- 1 - dat@prob[,idxNotDiff] #p(differentail genes)
   n <- nrow(dat@prob)
   fdr.hat <- sum(tau0 * (tau0 <= threshold)) / sum((tau0 <= threshold) ) 
   fndr.hat <- sum( (tau1) * (tau0 >= threshold))  / (n - sum(tau0 <= threshold) )
   fpr.hat <- sum(tau0 * (tau0 <= threshold)) / sum(tau0)
   fnr.hat <- sum(tau1 * (tau0> threshold) ) / sum(tau1)
return(list(FDR = fdr.hat, FNDR=fndr.hat, FPR = fpr.hat, FNR=fnr.hat) )
}


simDat <- function(mix){
    #simulate data for parametric bootstrap & compareMixModels
    #if weights are != 1, other parameters needs to be == 1!!! (--> otherwise, data cannot be unpacked reasonably!) 
     
    n= mix@num.obs #number of observations 

    k= mix@num.k
    myunif <- runif(n)
    z <- matrix(0, ncol=k, nrow=n)
    for (i in 1:k)
		if (i==1) z[,i] <- as.integer(sum(mix@p[1]) > myunif )
	    else z[,i] <- as.integer(sum(mix@p[1:(i-1)]) < myunif & (sum(mix@p[1:i]) >= myunif) )
    if (k==1) z <-matrix(1, ncol=k, nrow=n)
	#print (z)
    # z is a matrix that randomly assigns each observation (row) to a component (column) 
    # --> each row consists of 1x1 and (k-1)x0
    x = matrix(0,ncol=mix@num.k, nrow=n)

    if (mix@family == "gaussian"){
        if (sum(mix@dat[,4] == rep(1, n) ) == n) #dat[,4] only consists of ones
        tmpsd <-  sqrt(mix@component.var)
        else tmpsd <- mix@dat[,4]
        
        for (i in 1:k)
           x[,i] <- z[,i] * rnorm(n, mean=mix@t[i], sd= sqrt(tmpsd))
        res = apply(x,1,sum)  #a row consists of one number != 0 and (k-1) numbers == 0
    }
    else if (mix@family == "poisson"){
        Ei <- mix@dat[,3]
        for (i in 1:k)
           x[,i] <- z[,i] * rpois(n, mix@t[i]*Ei)
        res = apply(x,1,sum)  #a row consists of one number != 0 and (k-1) numbers == 0 
    }
    if (mix@family == "binomial"){
        tmpsize <- mix@dat[,3]
        for (i in 1:k)
           x[,i] <- z[,i] * rbinom(n, tmpsize, mix@t[i]) 
        res = apply(x,1,sum)  #a row consists of one number != 0 and (k-1) numbers == 0                 
    }
return(res)
}





summary.CAMAN.object <- function(object, ...){
	cl <- match.call()
	cat("Summary of a Computer Assisted Mixture Analysis: \n \n")
	cat("Data consists of", object@num.obs, "observations (rows). \n")
	cat("The Mixture Analysis identified", object@num.k, "component")
	if (length(object@num.k) >0) cat("s")
	cat(" of a", object@family, "distribution: \n \n")
	
	details <- matrix(0, nrow=object@num.k, ncol=2)
	descr_var <- ""
	if (object@family == "gaussian") {
		tmp <- "mean"
		if (object@is_metaAnalysis == 0) descr_var <- paste("component variance:", object@component.var, "\n")
	}
	else if (object@family == "poisson") tmp <- "lambda"
	else if (object@family == "binomial") tmp <- "prob"
	colnames(details) = c("p", tmp)
	rownames(details) = 1:object@num.k
	details[,1] <- object@p
	details[,2] <- object@t
	cat("DETAILS:\n")
	print(details)
	cat(descr_var)
	cat("\n \n")
	
	cat("Classification:\n")
    if (object@num.obs > 20 || object@num.k>8) cat("The classification matrix is too big to visualize here. \n type ",as.character(cl[2]),"@prob to watch the probability matrix \n or type ",as.character(cl[2]),"@classification to watch the \n class labeling (each row was assigned to its most likely component!)\n", sep="")
	else {
		cat("Classification matrix:\n")
		print(object@prob)
		cat("Class labeling:\n")
		print(object@classification)
	}
	cat("\n \nnumber of VEM-iterations done:",object@steps[1],"\n")        
	cat("alteration within the final VEM-iteration step:",object@finalacc[1],"\n")
	cat("number of EM-iterations done:",object@steps[2],"\n")        
	cat("alteration within the final EM-iteration step:",object@finalacc[2],"\n \n")
	cat("Log-Likelihood:",object@LL,"    ")
	cat("BIC:",object@BIC,"\n \n")
	
	# @otherParams = c(limit, numiter, acc, startk)
	
	cat("User-defined parameters:\n")
	cat("   max number of iterations:",object@otherParams[2],"\n")
	cat("   limit for combining components:",object@otherParams[1],"\n")
	cat("   threshold for converging:",object@otherParams[3],"\n")
	cat("   number of grid points (startk):",object@otherParams[4],"\n")
}


#some abbrevated commands
#mixboot <- mixalg.boot
#mixalg.Boot <- mixalg.boot
#mixpboot <- mixalg.paraBoot
#mix.anova <- anova.CAMAN.object
##########################################################################
#
#
# New bivariate functions
##########################################################################
#
bivariate.EM<-function(obs1,obs2,type,data = NULL, var1, var2, corr, lambda1, lambda2,p,numiter=5000,acc=1.e-7,class){

  ## avoid attach/detach but keep data argument optional
  if(!is.null(data)){
    cl <- match.call()
    varname <- function(x){sub("\\(\\)", "", deparse(x))}
    obs1_name <- varname(cl[2])
    obs2_name <- varname(cl[3])
    stopifnot(c(obs1_name, obs2_name) %in% names(data))
    obs1 <- getElement(data, obs1_name)
    obs2 <- getElement(data, obs2_name)

    var1_name <- varname(cl["var1"])
    var2_name <- varname(cl["var2"])
    if(!var1_name == "NULL" | !var2_name == "NULL"){
      stopifnot(c(var1_name, var2_name) %in% names(data)) 
      var1 <- getElement(data, var1_name)
      var2 <- getElement(data, var2_name)
    }

    corr_name <- varname(cl["corr"])
    if(!corr_name == "NULL"){
      stopifnot(corr_name %in% names(data)) 
      corr <- getElement(data, corr_name)
    }
  }
if(type=="bi"&& lambda1!=0 && lambda2 !=0 && p !=0){
##cat("### EM-algorithm for bivariate normally distributed data")


z1 <- function(a,n,l1,l2,pro, numiter,acc){.Call("ema_versh_st", as.vector(a), as.vector(n),as.vector(l1),as.vector(l2),as.vector(pro), as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er <- z1(obs1,obs2,lambda1,lambda2,p,numiter,acc)
len<-length(er)
l<-len/7
lam1<-er[1:l]
lam2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
var1<-er[((3*l)+1):(4*l)]
var2<-er[((4*l)+1):(5*l)]
corr<-er[((5*l)+1):(6*l)]
ll<-er[((6*l)+1)]
bic <- -2 * ll[1] + (3*length(prob)- 1) * log(length(obs1))

ERG<-matrix(data=c(lam1,lam2,prob,var1,var2,corr),nrow=l,ncol=6)
colnames(ERG) <- c("Lambda1","Lambda_2","Prob","Var1","Var2","Corr")
z2<-function(a,n,l1,l2,pro, numiter,acc){.Call("ema_ind_st", as.vector(a), as.vector(n),as.vector(l1),as.vector(l2),as.vector(pro), as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er2<-z2(obs1,obs2,lambda1,lambda2,p, numiter,acc)

x<-rep(1,629)
y<-rep(1,629)
a<-sqrt(qchisq(0.95,2))
t<-seq(0,6.28,0.01)
x<-array(rep(0,629*1*l),c(629,1,l))
y<-array(rep(0,629*1*l),c(629,1,l))
z<-array(rep(0,629*2*l),c(629,2,l))
for (i in 1:l){
x[, , i]<-lam1[i]+sqrt(var1[i])*a*cos(t)
y[,,i]<-lam2[i]+sqrt(var2[i])*a*cos(t+acos(corr[i]))
} 
for (i in 1:l){
z[, , i]<-c(x[,,i],y[,,i])
}
id<-er2
id<-id+1
a<-0
nn<-length(obs1)
mat<-matrix(data=c(obs1,obs2,id),nrow=nn,ncol=3)

if (class=="TRUE"){
res<-new("CAMAN.BIEM.object", RESULT=ERG ,BIC=bic,LL=ll[1],Mat=mat, Z=z,cl=er2)}
if (class=="FALSE"){
res<-new("CAMAN.BIEM.object", RESULT=ERG ,BIC=bic,LL=ll[1],Mat=mat, Z=z)}

return(res)
}

if(type=="meta" &&  lambda1!=0 && lambda2!=0 && p!=0){

z3 <- function(a,n,v1,v2,l1,l2,pro, numiter,acc){
  .Call("ema_meta_st", as.vector(a), as.vector(n),as.vector(v1),as.vector(v2),
        as.vector(l1),as.vector(l2),as.vector(p),as.integer(numiter),as.double(acc), 
        PACKAGE = "CAMAN")}
er<-z3(obs1,obs2,var1,var2,lambda1,lambda2,p,numiter,acc)
len<-length(er)
l<-len/5
lam1<-er[1:l]
lam2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
ll<-er[((3*l)+1)]
max_grad<-er[(4*l+1):len]
bic <- -2 * ll[1] + (3*length(prob)- 1) * log(length(obs1))

ERG1<-matrix(data=c(lam1,lam2,prob),nrow=l,ncol=3)
colnames(ERG1) <- c("lambda1","lambda2","p")
#colnames(mat) <- c("Lambda_1","Lambda_2","Prob","LL","max_grad")
z4<-function(a,n,v1,v2,l1,l2,p, numiter,acc){.Call("ema_ind_meta_st", as.vector(a), as.vector(n),as.vector(v1),as.vector(v2),as.vector(l1),as.vector(l2),as.vector(p),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z4(obs1,obs2,var1,var2,lambda1,lambda2,p, numiter,acc)



if (class=="TRUE"){
res<-new("CAMAN.BIEM.object", RESULT=ERG1 ,BIC=bic,LL=ll[1], cl=er)}
if (class=="FALSE"){
res<-new("CAMAN.BIEM.object", RESULT=ERG1 ,BIC=bic,LL=ll[1])}

return(res)

}
} 
bivariate.VEM <-function(obs1,obs2,type,data = NULL, var1, var2, lambda1, lambda2,p, startk, numiter=5000,acc=1.e-7){

  ## avoid attach/detach but keep data argument optional
  if(!is.null(data)){
    cl <- match.call()
    varname <- function(x){sub("\\(\\)", "", deparse(x))}
    obs1_name <- varname(cl[2])
    obs2_name <- varname(cl[3])
    stopifnot(c(obs1_name, obs2_name) %in% names(data))
    obs1 <- getElement(data, obs1_name)
    obs2 <- getElement(data, obs2_name)
    
    var1_name <- varname(cl["var1"])
    var2_name <- varname(cl["var2"])
    if(!var1_name == "NULL" | !var2_name == "NULL"){
      stopifnot(c(var1_name, var2_name) %in% names(data)) 
      var1 <- getElement(data, var1_name)
      var2 <- getElement(data, var2_name)
    }
  }
  
  if(type=="uni"){
z5<-function(a, startk, numiter,acc){.Call("vem_uni", as.vector(a),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z5(obs1, startk,numiter,acc)
len<-length(er)
l<-len/3
lam<-er[1:l]
prob<-er[(l+1):(2*l)]
ll<-er[((2*l)+1)]
bic <- -2 * ll[1] + (3*2- 1) * log(length(obs1))
mat<-matrix(data=c(lam,prob),nrow=l,ncol=2)
colnames(mat) <- c("lambda","mixing Prob")
res<-new("CAMAN.BIVEM.object", RESULT_uni=mat ,BIC=bic,LL=ll[1])
return(res)

}
if(type=="bi"){

z6<-function(a,n, startk, numiter,acc){.Call("vem_bi_sh", as.vector(a), as.vector(n),as.integer(startk), as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z6(obs1,obs2, startk,numiter,acc)
len<-length(er)
l<-len/4
lam1<-er[1:l]
lam2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
ll<-er[((3*l)+1)]
bic <- -2 * ll[1] + (3*length(prob)- 1) * log(length(obs1))
mat<-matrix(data=c(lam1,lam2,prob),nrow=l,ncol=3)
#colnames(mat) <- c("Lambda_1","Lambda_2","Prob")

res<-new("CAMAN.BIVEM.object", RESULT=mat ,BIC=bic,LL=ll[1])
return(res)
#res<-list("Vem for bivariate data","Lambda_1"=lam1, #"Lambda_2"=lam2,"Prob"=prob)
#res<-list("VEM algorithm for bivariate data",mat)
#print(res)
#cat("BIC : ", bic,"\n")
#cat("Log-Likelihood: ", ll[1],"\n")
} 


if(type=="meta"){
z7<-function(a,n,v1,v2, startk,numiter,acc){.Call("vem_versh_meta_sh", as.vector(a), as.vector(n),as.vector(v1),as.vector(v2),as.integer(startk), as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z7(obs1,obs2,var1,var2, startk,numiter,acc)
len<-length(er)
l<-len/4
lam1<-er[1:l]
lam2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
ll<-er[((3*l)+1) :(4*l)]
bic <- -2 * ll[1] + (3*3- 1) * log(length(obs1))
Mat<-matrix(data=c(prob,lam1,lam2),nrow=l,ncol=3)
#colnames(mat) <- c("p","lambda1","lambda1")
 
#res<-list("VEM algorithm for diagnostic meta analysis", mat)
res<-new("CAMAN.BIVEM.object", RESULT_meta=Mat,BIC=bic,LL=ll[1])
return(res)
#res<-list("VEM algorithm for diagnostic meta analysis", "lambda_1"=lam1, "lambda_2"=lam2,"p"=prob)
#print(res)
#cat("BIC : ", bic,"\n")
#cat("Log-Likelihood: ", ll[1],"\n")

}
}

vem_grad<-function(obs1,obs2,type,data = NULL,var1, var2,lambda1, lambda2,p, startk,numiter=5000,acc=1.e-7){

  ## avoid attach/detach but keep data argument optional
  if(!is.null(data)){
    cl <- match.call()
    varname <- function(x){sub("\\(\\)", "", deparse(x))}
    obs1_name <- varname(cl[2])
    obs2_name <- varname(cl[3])
    stopifnot(c(obs1_name, obs2_name) %in% names(data))
    obs1 <- getElement(data, obs1_name)
    obs2 <- getElement(data, obs2_name)
    
    var1_name <- varname(cl["var1"])
    var2_name <- varname(cl["var2"])
    if(!var1_name == "NULL" | !var2_name == "NULL"){
      stopifnot(c(var1_name, var2_name) %in% names(data)) 
      var1 <- getElement(data, var1_name)
      var2 <- getElement(data, var2_name)
    }
  }
  
z14<-function(a,n, startk, numiter,acc){.Call("vem_bi_grad", as.vector(a), as.vector(n),as.integer(startk), as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z14(obs1,obs2, startk,numiter,acc)
len<-length(er)
l<-len/4
lam1<-er[1:l]
lam2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
grad<-er[((3*l)+1):(4*l)]
#bic <- -2 * ll[1] + (3*5- 1) * log(length(obs1))
mat<-matrix(data=c(lam1,lam2,prob,grad),nrow=l,ncol=4)
colnames(mat) <- c("Lambda_1","Lambda_2","Prob","Grad")
#res<-list("Vem for bivariate data","Lambda_1"=lam1, "Lambda_2"=lam2,"Prob"=prob, "Grad"=grad)
res<-list("VEM algorithm for bivariate data",mat)
print(res)
}


bivariate.mixalg<-function(obs1,obs2,type,data = NULL,var1, var2, corr, lambda1, lambda2,p,startk, numiter=5000,acc=1.e-7,class){

  ## avoid attach/detach but keep data argument optional
  if(!is.null(data)){
    cl <- match.call()
    varname <- function(x){sub("\\(\\)", "", deparse(x))}
    obs1_name <- varname(cl[2])
    obs2_name <- varname(cl[3])
    stopifnot(c(obs1_name, obs2_name) %in% names(data))
    obs1 <- getElement(data, obs1_name)
    obs2 <- getElement(data, obs2_name)
    
    var1_name <- varname(cl["var1"])
    var2_name <- varname(cl["var2"])
    if(!var1_name == "NULL" | !var2_name == "NULL"){
      stopifnot(c(var1_name, var2_name) %in% names(data)) 
      var1 <- getElement(data, var1_name)
      var2 <- getElement(data, var2_name)
    }
    corr_name <- varname(cl["corr"])
    if(!corr_name == "NULL"){
      stopifnot(corr_name %in% names(data)) 
      corr <- getElement(data, corr_name)
    }
  }
  
  if(type=="uni"){
z8<-function(a, startk, numiter,acc){.Call("ema_uni", as.double(a),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z8(obs1, startk, numiter,acc)
len<-length(er)
l<-len/3
lam<-er[1:l]
prob<-er[(l+1):(2*l)]
var<-er[(2*l+1):len]
matu<-matrix(data=c(lam,prob,var),nrow=l,ncol=3)
#colnames(mat) <- c("Lambda","Prob","Var")
z9<-function(a,startk, numiter,acc){.Call("ema_ind_uni", as.double(a),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z9(obs1,startk, numiter,acc)

if (class=="TRUE"){
res<-new("CAMAN.BIMIXALG.object", RESULT_uni=matu ,BIC=bic,LL=ll[1],cl=er)}
if (class=="FALSE"){
res<-new("CAMAN.BIMIXALG.object", RESULT_uni=matu ,BIC=bic,LL=ll[1])}

return(res)


}

if(type=="bi"&& lambda1==0 && lambda2 ==0 && p==0){

z10<-function(a,n,startk, numiter,acc){.Call("ema_versh_sh", as.double(a), as.double(n),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z10(obs1,obs2,startk, numiter,acc)
len<-length(er)
l<-len/7
lambda1<-er[1:l]
lambda2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
var1<-er[((3*l)+1):(4*l)]
var2<-er[((4*l)+1):(5*l)]
corr<-er[((5*l)+1):(6*l)]
ll<-er[((6*l)+1)]
bic <- -2 * ll[1] + (3*length(prob)- 1) * log(length(obs1))
ERG<-matrix(data=c(lambda1,lambda2,prob,var1,var2,corr),nrow=l,ncol=6)
#colnames(ERG) <- c("Lambda1","Lambda_2","Prob","Var1","Var2","Corr")
z11<-function(a,n,startk, numiter,acc){.Call("ema_ind_sh", as.double(a), as.double(n),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er2<-z11(obs1,obs2,startk, numiter,acc)
#par(mfrow=c(2,1))
#plot(obs1,obs2, xlab = "x1", ylab = "x2",pch=19,col="blue",cex=0.4,main=" rs12363681");
x<-rep(1,629)
y<-rep(1,629)
a<-sqrt(qchisq(0.95,2))
t<-seq(0,6.28,0.01)
x<-array(rep(0,629*1*l),c(629,1,l))
y<-array(rep(0,629*1*l),c(629,1,l))
z<-array(rep(0,629*2*l),c(629,2,l))
for (i in 1:l){
x[, , i]<-lambda1[i]+sqrt(var1[i])*a*cos(t)
y[,,i]<-lambda2[i]+sqrt(var2[i])*a*cos(t+acos(corr[i]))
}
for (i in 1:l){
z[, , i]<-c(x[,,i],y[,,i])
}
id<-er2
id<-id+1
a<-0
nn<-length(obs1)
mat<-matrix(data=c(obs1,obs2,id),nrow=nn,ncol=3)

if (class=="TRUE"){
res<-new("CAMAN.BIMIXALG.object", RESULT=ERG ,BIC=bic,LL=ll[1],Mat=mat, Z=z,cl=er2)}
if (class=="FALSE"){
res<-new("CAMAN.BIMIXALG.object", RESULT=ERG ,BIC=bic,LL=ll[1],Mat=mat, Z=z)}

return(res)

}
if(type=="meta" &&  lambda1==0 && lambda2==0 && p==0){

z12<-function(a,n,v1,v2, startk, numiter,acc){.Call("ema_meta_sh", as.double(a), as.double(n),as.vector(v1),as.vector(v2),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z12(obs1,obs2,var1,var2, startk, numiter,acc)
len<-length(er)
l<-len/5
lambda1<-er[1:l]
lambda2<-er[(l+1):(2*l)]
prob<-er[((2*l)+1):(3*l)]
ll<-er[(3*l+1):(4*l)]
max_grad<-er[(4*l+1):len]
bic <- -2 * ll[1] + (3*length(prob)- 1) * log(length(obs1))

ERG1<-matrix(data=c(lambda1,lambda2,prob),nrow=l,ncol=3)
#colnames(ERG1) <- c("Lambda_1","Lambda_2","Prob")
z13<-function(a,n,v1,v2,startk, numiter,acc){.Call("ema_ind_meta_sh", as.double(a), as.double(n),as.vector(v1),as.vector(v2),as.integer(startk),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
er<-z13(obs1,obs2,var1,var2,startk, numiter,acc)

if (class=="TRUE"){
res<-new("CAMAN.BIMIXALG.object", RESULT_meta=ERG1 ,BIC=bic,LL=ll[1],cl=er)}
if (class=="FALSE"){
res<-new("CAMAN.BIMIXALG.object", RESULT_meta=ERG1 ,BIC=bic,LL=ll[1],)}
return(res)

}
}



CAMANboot<-function(obs1,obs2,var1,var2,lambda11,lambda12,prob1,lambda21,lambda22,prob2,rep,data,numiter=10000,acc=1.e-7){

    cl <- match.call()
    varname <- function(x){sub("\\(\\)", "", deparse(x))}
    obs1_name <- varname(cl[2])
    obs2_name <- varname(cl[3])
    stopifnot(c(obs1_name, obs2_name) %in% names(data))
    obs1 <- getElement(data, obs1_name)
    obs2 <- getElement(data, obs2_name)
    
    var1_name <- varname(cl["var1"])
    var2_name <- varname(cl["var2"])
    if(!var1_name == "NULL" | !var2_name == "NULL"){
      stopifnot(c(var1_name, var2_name) %in% names(data)) 
      var1 <- getElement(data, var1_name)
      var2 <- getElement(data, var2_name)
    }
  
  a<-matrix(nrow=(length(data[,1])),ncol=2)
k_1<-matrix(nrow=rep,ncol=(length(lambda11)*3+2))
k_2<-matrix(nrow=rep,ncol=(length(lambda21)*3+2))

#print("###Bootstrap fuer Metadaten mit Startwerten")
fun<-function(a,n,v1,v2,l1,l2,pro,numiter,acc){.Call("ema_meta_st", as.vector(a), as.vector(n),as.vector(v1),as.vector(v2),as.vector(l1),as.vector(l2),as.vector(pro),as.integer(numiter), as.double(acc), PACKAGE = "CAMAN")}
j<-0
repeat{
j<-j+1
for(i in 1:(length(data[,1]))){
z<-runif(1)
size<-length(prob1)
if (size>1){
cdf<-rep(1:(size-1))
for(m in 1:(size-1)){
cdf[1]<-prob1[1]
cdf[m+1]<-(cdf[m]+prob1[m+1])}
if (z<cdf[1]){
a[i,]<-rmvnorm(n = 1, mean=c(lambda11[1],lambda12[1]),sigma <- matrix(c(var1[i],0,0,var2[i]), ncol=2))
}
for(ll in 2:(size)){
if (z>cdf[ll-1]&&z<cdf[ll]){
a[i,]<-rmvnorm(n =1, mean=c(lambda11[ll],lambda12[ll]),sigma <- matrix(c(var1[i],0,0,var2[i]), ncol=2))

if (ll<(size-1)){
ll=ll+1}
else break
rm(ll)
}
}
if (z>=cdf[size-1]){
a[i,]<-rmvnorm(n =1, mean=c(lambda11[size],lambda12[size]),sigma <- matrix(c(var1[i],0,0,var2[i]), ncol=2))
}
}
if (size==1){
a[i,]<-rmvnorm(n = 1, mean=c(lambda11[1],lambda12[1]),sigma <- matrix(c(var1[i],0,0,var2[i]), ncol=2))
}
}
er1<-fun(a[,1],a[,2],var1,var2,lambda11,lambda12,prob1,numiter,acc)
er2<-fun(a[,1],a[,2],var1,var2,lambda21,lambda22,prob2,numiter,acc)

len1<-length(er1)
l1<-len1/5
lambda1_1<-er1[1:l1]
lambda1_2<-er1[(l1+1):(2*l1)]
prob_1<-er1[((2*l1)+1):(3*l1)]
ll_1<-er1[(3*l1+1)]
max_grad_1<-er1[(4*l1+1)]
len2<-length(er2)
l2<-len2/5
lambda2_1<-er2[1:l2]
lambda2_2<-er2[(l2+1):(2*l2)]
prob_2<-er2[((2*l2)+1):(3*l2)]
ll_2<-er2[(3*l2+1)]
max_grad_2<-er2[(4*l2+1)]
#k_1<-matrix(nrow=3,ncol=((l*4)-1))
k_1[j,]<-matrix(data=c(lambda1_1,lambda1_2,prob_1,ll_1,max_grad_1),nrow=1,ncol=((length(lambda1_2)*3)+2))
k_2[j,]<-matrix(data=c(lambda2_1,lambda2_2,prob_2,ll_2,max_grad_2),nrow=1,ncol=((length(lambda2_2)*3)+2))
k_1<-round(x=k_1,digits=3)
k_2<-round(x=k_2,digits=3)
m1<-length(lambda1_2)*3+1
m2<-length(lambda2_2)*3+1
ii<-1:length(lambda1_1)
aa<-paste("lam1",ii,sep="")
bb<-paste("lam2",ii,sep="")
cc<-paste("prob",ii,sep="")
colnames(k_1)<-c(aa,bb,cc,"LL_1","max_grad")
i<-1:length(lambda2_1)
aaa<-paste("lam1",i,sep="")
bbb<-paste("lam2",i,sep="")
ccc<-paste("prob",i,sep="")
colnames(k_2)<-c(aaa,bbb,ccc,"LL_2","max_grad")
#colnames(k_1)<-c("lambda_1","lambda_2","prob","ll")
if(j==rep) break
}

var_1<-rep(1,TRUE)
var_2<-rep(1,TRUE)
corr_1<-rep(1,TRUE)
covv_1<-rep(1,TRUE)
corr_2<-rep(1,TRUE)
covv_2<-rep(1,TRUE)
S1 <- array(rep(0, 2 * 2 * l1), c(2, 2, l1))
S2 <- array(rep(0, 2 * 2 * l2), c(2, 2, l2))
for(i in 1:(l1*2)){
var_1[i]<-var(k_1[,i], use = "complete")

}

m<-1

for(i in 1:l1 ){
corr_1[i]<-cor(k_1[,m],k_1[,(m+1)], use = "complete")
covv_1[i]<-(sqrt(var_1[m])*sqrt(var_1[(m+1)]))*corr_1[i]

m=m+2}


m<-1
for(i in 1:l1 ){
S1[, , i]<-rbind(var_1[m],covv_1[i],covv_1[i],var_1[(m+1)])
m=m+2
}

m<-1
for(i in 1:(l2*2)){

var_2[i]<-var(k_2[,i], use = "complete")}

for(i in 1:l2){
corr_2[i]<-cor(k_2[,m],k_2[,(m+1)], use = "complete")
covv_2[i]<-(sqrt(var_2[m])*sqrt(var_2[m+1]))*corr_2[i]
m=m+2}
 m<-1
for(i in 1:l2 ){
S2[, , i]<-rbind(var_2[m],covv_2[i],covv_2[i],var_2[(m+1)])
m=m+2
}

for(i in rep){
differenz<-(k_1[,m1]-k_2[,m2])}
llh<-(-2)*differenz
s_ll<-sort(llh)


res<-new("CAMAN.BOOT.object", H0=k_1, S1=S1, H1=k_2,S2=S2, LL=s_ll, Q95=quantile(s_ll,0.95), Q975=quantile(s_ll,0.975), Q99=quantile(s_ll,0.99))
return(res)
} 
#some abbrevated commands
mixboot <- mixalg.boot
mixalg.Boot <- mixalg.boot
mixpboot <- mixalg.paraBoot
mix.anova <- anova.CAMAN.object
