####Function to run the sampler#####
##SPECIFYING PRIORS######
#########################
##if priors = NULL => uses randomly generated priors
##else: priors is a list of following objects:
        ##MuBeta:= prior mean for betas & intercepts; 
        ##SigmaBeta:= prior variance for betas & intercept;
        ##MuAlpha:= prior mean for alpha; 
        ##SigmaAlpha:= prior variance for alpha;
        ##MuZ; VarZ;
        ##PriorA; PriorB 
##TUNING PARAMETERS#######
##########################
##if tune = NULL => uses auto tuning
##else: tune is a list of following objects:
       ##tuneAlpha = 0.9
       ##tuneBeta = vec, PP
       ##tuneInt = vec, len = 1
       ##tuneZ =  list( vec(len = nn[x]])) length of list = KK
##############################################################          
#library(MASS)
HLSMfixedEF= function(Y,edgeCov = NULL, receiverCov = NULL,senderCov =NULL,FullX = NULL, initialVals = NULL, priors = NULL, tune = NULL,
        tuneIn = TRUE, TT = NULL,dd, niter,intervention)
{

  #X and Y are provided as list. 
    if(class(Y) != 'list'){
	if(dim(Y)[2] != 4){stop('Invalid data structure type')} }

    if(class(Y) == 'list' & class(Y[[1]]) != 'matrix' & class(Y[[1]]) != 'data.frame'){stop('Invalid data structure type')}	
    if(class(Y) == 'list'){ 
        KK = length(Y)
	check = 1
	for(cc in 1:KK){
	    check =  check *(dim(Y[[cc]])[1] == dim(Y[[cc]])[2])
	}
	if(check == 1){	
	    nn =sapply(1:KK,function(x) nrow(Y[[x]])) 
	    nodenames = lapply(1:KK,function(x) dimnames(Y[[x]])[[1]])
		}

	if(check == 0 & dim(Y[[1]])[2] == 4){
		nn = sapply(1:KK, function(x)length(unique(c(Y[[x]]$Receiver,Y[[x]]$Sender))))
		nodenames = lapply(1:KK, function(x) unique(c(Y[[x]]$Receiver,Y[[x]]$Sender)))
	}	}

    if(class(Y) != 'list'){
	if(dim(Y)[2] == 4){
		nid = unique(Y$id)
		KK = length(nid)
		nn = rep(0,KK)
		df.list = list()
		nodenames = list()
		for(k in 1:KK){
			df.sm = Y[which(Y$id == nid[k],),]
			nn[k] = length(unique(c(df.sm$Receiver,df.sm$Sender)))
			nodenames[[k]] = unique(c(df.sm$Receiver, df.sm$Sender))
			df.list[[k]] = array(0, dim = c(nn[k],nn[k]))
			dimnames(df.list[[k]])[[1]] = dimnames(df.list[[k]])[[2]] = nodenames[[k]]
			for(i in 1:dim(df.sm)[1]){
				df.list[[k]][df.sm$Sender[i],df.sm$Receiver[i]] = df.sm$Outcome[i]  #assume undirected graph and missing items are zeros
			}
		}
		Y = df.list 
	}}


##prepare covariates#####
#########################
	noCOV = FALSE
	if(!is.null(FullX) & !is.null(edgeCov) &!is.null(receiverCov) & !is.null(senderCov))(stop('FullX cannot be used when nodal or edge covariates are provided'))

	if(is.null(FullX) & is.null(edgeCov) & is.null(receiverCov) & is.null(senderCov)){
		X = lapply(1:KK,function(x) array(0, dim = c(nn[x],nn[x],1)))
		noCOV = TRUE
	}

	if(is.null(FullX)){
	if(!is.null(edgeCov) | !is.null(senderCov)| !is.null(receiverCov)){
	  if(!is.null(edgeCov)){
		if(class(edgeCov) != 'data.frame'){
			stop('edgeCov must be of class data.frame')}
		X1 = getEdgeCov(edgeCov, nn,nodenames)
}else(X1 =NULL)
  	  if(!is.null(senderCov)){
		if(class(senderCov) != 'data.frame'){
			stop('senderCov must be of class data.frame')}
		X2 = getSenderCov(senderCov, nn,nodenames)
}else(X2 = NULL)


	  if(!is.null(receiverCov)){
		if(class(receiverCov) != 'data.frame'){
			stop('receiverCov must be of class data.frame')}
		X3 = getReceiverCov(receiverCov, nn,nodenames)
}else(X3 = NULL)	

	X = lapply(1:KK, function(x){if(!is.null(X1)&!is.null(X2)&!is.null(X3)){
		ncov = dim(X1[[x]])[3]+dim(X2[[x]])[3]+dim(X3[[x]])[3];
		df = array(0, dim = c(nn[x],nn[x],ncov));
		df[,,1:dim(X1[[x]])[3]] = X1[[x]];
		df[,,(dim(X1[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X2[[x]])[3])] = X2[[x]];
		df[,,(dim(X1[[x]])[3]+dim(X2[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X2[[x]])[3]+dim(X3[[x]])[3])] = X3[[x]] };
		if(!is.null(X1)&!is.null(X2) & is.null(X3)){
			ncov = dim(X1[[x]])[3]+dim(X2[[x]])[3];
			df = array(0, dim = c(nn[x],nn[x],ncov));
			df[,,1:dim(X1[[x]])[3]] = X1[[x]];
			df[,,(dim(X1[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X2[[x]])[3])] = X2[[x]]};
		if(!is.null(X1)&!is.null(X3)&is.null(X2)){
			ncov = dim(X1[[x]])[3]+dim(X3[[x]])[3];
			df = array(0, dim = c(nn[x],nn[x],ncov));
			df[,,1:dim(X1[[x]])[3]] = X1[[x]];
			df[,,(dim(X1[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X3[[x]])[3])] = X3[[x]]};
	if(!is.null(X2)&!is.null(X3)&is.null(X1)){
			ncov = dim(X2[[x]])[3]+dim(X3[[x]])[3];
			df = array(0, dim = c(nn[x],nn[x],ncov));
			df[,,1:dim(X2[[x]])[3]] = X2[[x]];
			df[,,(dim(X2[[x]])[3]+1):(dim(X2[[x]])[3]+dim(X3[[x]])[3])] = X3[[x]]};
	if(!is.null(X1)& is.null(X2)& is.null(X3)){
			df = X1[[x]] };
	if(is.null(X1)& !is.null(X2)& is.null(X3)){
			df = X2[[x]] };
	if(is.null(X1)& is.null(X2)& !is.null(X3)){
			df = X3[[x]] };
	return(df) } )
}
}
	if(!is.null(FullX)) X = FullX


#    nn = sapply(1:length(X),function(x) nrow(X[[x]]))
#    KK = length(X)
    PP = dim(X[[1]])[3]
    XX = unlist(X)
    YY = unlist(Y)
    YY[which(is.na(YY))] = 0
    XX[which(is.na(XX))] = 0
    #Priors

    if(is.null(priors)){
	MuBeta= rep(0,(PP+1)) 
	VarBeta = rep(1,(PP+1)) 
        MuAlpha=0 
        VarAlpha = 1 
        MuZ = c(0,0)
        VarZ = c(20,20)
        PriorA = 100
        PriorB = 150
     }else{
	if(class(priors) != 'list')(stop("priors must be of class list, if not NULL"))
	MuBeta = priors$MuBeta
	VarBeta = priors$VarBeta
	MuAlpha = priors$MuAlpha
	VarAlpha = priors$VarAlpha
	MuZ = priors$MuZ
	VarZ = priors$VarZ
	PriorA = priors$PriorA
	PriorB = priors$PriorB
  }
##starting values
    if(is.null(initialVals)){
	Z0 = list()
        for(i in 1:KK){  
            ZZ = t(replicate(nn[i],rnorm(dd,0,1)))
            ZZ[1,]=c(1,0)
            ZZ[2,2]=0
            if(ZZ[2,1] < ZZ[1,1]){
                ZZ[2,1] = -1*(ZZ[2,1]-ZZ[1,1])+1}
            ZZ[3,2] = abs(ZZ[3,2])
	    Z0[[i]] = ZZ		 
    }
        Z0 = unlist(Z0)
        beta0 = rnorm(PP,0,1)
        intercept0  = rnorm(1, 0,1)
        if(intervention == 1){    alpha0=rnorm(1, 0, 1) }
        print("Starting Values Set")
    }else{
	if(class(initialVals) != 'list')(stop("initialVals must be of class list, if not NULL"))
	Z0 = initialVals$ZZ
	beta0 = initialVals$beta
	intercept0 = initialVals$intercept
	if(intervention == 1){ alpha0 = initialVals$alpha}
	}

    if(intervention == 0){
        alpha0 = 0
        TT = rep(0, KK)
    }

###tuning parameters#####
    if(is.null(tune)){
            a.number = 5
            tuneAlpha = 0.9
            tuneBeta = rep(1,PP)
            tuneInt = 0.2
            tuneZ =  lapply(1:KK,function(x) rep(1.2,nn[x]))          
            } else{
         	if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
                 a.number = 1
                 tuneAlpha = tune$tuneAlpha
                 tuneBeta = tune$tuneBeta
                 tuneInt = tune$tuneInt
                 tuneZ = tune$tuneZ
          }       
  
###Tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
    while(do.again ==1){
        print('Tuning the Sampler')
        for(counter in 1:a.number)
{ 
            rslt = MCMCfixedEF(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,TT = TT,
		beta = beta0 ,intercept = intercept0,alpha = alpha0,MuAlpha = MuAlpha,SigmaAlpha = VarAlpha,
		MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt,
		tuneAlpha = tuneAlpha,tuneZAll = unlist(tuneZ),niter = 200,PriorA = PriorA, PriorB = PriorB, 
		intervention = intervention)

            tuneAlpha = adjust.my.tune(tuneAlpha, rslt$acc$alpha,1)
     	    tuneZ = lapply(1:KK,function(x)adjust.my.tune(tuneZ[[x]], rslt$acc$Z[[x]], 2))
            tuneBeta = adjust.my.tune(tuneBeta, rslt$acc$beta,1)
            tuneInt =  adjust.my.tune(tuneInt,rslt$acc$intercept,1)
            print(paste('TuneDone = ',tuneX))
            tuneX = tuneX+1
    }
    extreme = lapply(1:KK,function(x)which.suck(rslt$acc$Z[[x]],2))
    do.again = max(sapply(extreme, length)) > 5
 
}
    print("Tuning is finished")  
}

    rslt = MCMCfixedEF(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,TT = TT,
		beta = beta0 ,intercept = intercept0,alpha = alpha0,MuAlpha = MuAlpha,SigmaAlpha = VarAlpha,
		MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt, 
		tuneAlpha = tuneAlpha,tuneZAll = unlist(tuneZ),niter = niter,PriorA = PriorA, PriorB = PriorB, 
		intervention = intervention)

    rslt$call = match.call()
    if(noCOV == TRUE & intervention == 0){
		rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt)	
		rslt$draws$Beta = NA
		rslt$draws$Alpha = NA
    }
    if(noCOV == TRUE & intervention == 1){
	    rslt$tune = list(tuneAlpha = tuneAlpha, tuneZ = tuneZ,tuneInt = tuneInt)
	    rslt$draws$Beta = NA	
	}
    if(noCOV == FALSE & intervention == 0){
	    rslt$tune = list(tuneBeta = tuneBeta, tuneZ = tuneZ,tuneInt = tuneInt)
	    rslt$draws$Alpha = NA
}
    if(noCOV == FALSE & intervention == 1){
	    rslt$tune = list(tuneBeta = tuneBeta,tuneAlpha=tuneAlpha,tuneZ = tuneZ,tuneInt = tuneInt)
}	
    class(rslt) = 'HLSM'
    rslt
}


    



