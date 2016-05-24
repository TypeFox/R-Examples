# TODO: Add comment
# 
# Author: Giorgio Spedicato
###############################################################################




cwrEm<-function(X,Y, nc, max_iter=1000, thresh=0.01, cov_typeX="full", cov_typeY="full",clamp_weights=FALSE, create_init_params=TRUE, cwrStart=NULL, cov_priorX=NULL,cov_priorY=NULL, verbose=TRUE, regress=TRUE, clamp_covX=FALSE, clamp_covY=FALSE) 
{
	#require("matlab")
	pmt<-proc.time();  #start process time monitor
	X=t(as.matrix(X)); #transpose matrix
	Y=t(as.matrix(Y)); #transpose matrix
	cwrOutput=list();  #creates output list
	dimensionX=size(X); #dim(X); 
	dimensionY=size(Y); #size(Y)
	nx=dimensionX[1];  #vars vettore X
	ny=dimensionY[1];  #vars vettore Y
	N=dimensionX[2];   #dimensione di X
	
	if (dimensionX[2] !=dimensionY[2]) stop("Unequal length of X and Y observations") #controls
	if (dimensionX[2] < dimensionX[1]) warning("More variables than observations in X")
	if (dimensionY[2] < dimensionY[1]) warning("More variables than observations in Y")
	if (nc > N) stop("more clusters than observations")
	
	if((create_init_params==FALSE) & (is.null(cwrStart))) stop("Needed suitable starting values of cwr"); #controls
	
	if (nc==1) { #works only if there is one cluster
		w=1/N;   #weights 
		WYbig=Y*w; #1/N*Y
		WYY=WYbig%*%t(Y); #E(Y^2)
		WY=as.matrix(apply(WYbig, 1, sum)); #is EY 
		WYTY = as.matrix(sum(diag(t(WYbig)%*%Y)));   #FIX sum(diag(t(WYbig)%*%Y));    may be a part of regression estimation equation
		cwrOutput$priorC=as.matrix(1);     #the prior probability is 1 (only one grou<e8>p
		cwrOutput$SigmaX=numeric(0);
		if (!regress) {
			cwrOutput$weightsY=numeric(0); #weightsY=beta of regression
			momentsY=.mixGaussMstep(1, WY, WYY, WYTY, cov_type=cov_typeY, cov_prior=cov_priorY); #fixes regression
			cwrOutput$muY=momentsY$mu;cwrOutput$SigmaY=momentsY$SigmaY;
			#qua ci vorrebbe l'assert
		} else {
			WXbig = X*w; #x ponderato
			WXX = WXbig%*%t(X); # #matrix E(XX)
			WX = as.matrix(apply(WXbig, 1, sum)); #matrix E(X)
			WXTX = as.matrix(sum(diag(t(WXbig)%*%X))); #fixed
			WXY = WXbig%*%t(Y); #these are regression equations;
			#calls weightsY (beta), muY (intercept), sig(sigma^2 epsilon) estimatos
			clgMoments=.clgMstep(1, WY, WYY, WYTY, WX, WXX, WXY, 
					cov_type= cov_typeY, cov_prior=cov_priorY); #estimaes muY and weights
			cwrOutput$muY=clgMoments$mu;cwrOutput$SigmaY=clgMoments$SigmaY; cwrOutput$weightsY=clgMoments$B; #save ouptuts
		}
		if(clamp_covY==TRUE) cwrOutput$SigmaY=cwrStart$SigmaY       #modified: recheck
		if(clamp_weights==TRUE) cwrOutput$weightsY=cwrStart$weightsY
		class(cwrOutput)="cwrObj" #assign the class
		return(cwrOutput)   #return fitted cwrObj and program ends
	}  else {  #nc>=2
		if(create_init_params){
			#initials parameters are calculated separately for X and Y (assuming therefore joint independence)
			initialsX = .mixgauss_init(nc, X, cov_typeX); #calculates initials fo X
			cwrOutput$muX=initialsX$mu;  #use kmeans
			cwrOutput$SigmaX=initialsX$Sigma;
			initialsY = .mixgauss_init(nc, Y, cov_typeY);  #calculate initials for Y
			cwrOutput$muY=initialsY$mu;
			cwrOutput$SigmaY=initialsY$Sigma;
			cwrOutput$weightsY = array(0, dim=c(ny, nx, nc));   #regression weights initial = 0 (beta)
			#regression weight defines P(y|x, Q) they lies in R^(dimY+dimX+1)
			cwrOutput$priorC=as.matrix(rep(1,nc)/nc)  #equal weights of group
		} else { #if are given it takes estimate parameters from cwrStart
			cwrOutput$muX=cwrStart$muX;
			cwrOutput$muY=cwrStart$muY;
			cwrOutput$SigmaX=cwrStart$SigmaX;
			cwrOutput$SigmaY=cwrStart$SigmaY;
			cwrOutput$weightsY = cwrStart$weightsY;
			cwrOutput$priorC = cwrStart$priorC;
		} 
	}
	if(clamp_covY==TRUE) cwrOutput$SigmaY = cwrStart$SigmaY   #they might have given  SigmaX, SigmaY and beta (to be checked)
	if(clamp_covX==TRUE)  cwrOutput$SigmaX = cwrStart$SigmaX #(to be checked)
	if(clamp_weights==TRUE)  cwrOutput$weightsY = cwrStart$weightsY #(to be checked)
	
	previous_logLik = -Inf;   #initialize estimation
	num_iter = 1;
	converged = FALSE;
	
	while ((num_iter <= max_iter) & (converged==FALSE)) #estimation 
	{
#   E step
		logLikList=.cwr_prob(cwrOutput, X, Y);  #estimate a list of probabiliies defined by (actual) expected statistics
		logLik=rowSums(log(logLikList$likXandY)); #this is the total likelihood of element points (sum of total likelihood)
		posteriors=logLikList$post #estimate posterior probabilies of pertaining to cluster i
		w = as.matrix(rowSums(posteriors)) #recalculates priors
		
		WYY = zeros(ny, ny, nc); #start conditional moment matrices
		WY = zeros(ny, nc); #
		WYTY = zeros(nc,1); #
		
		WXX = zeros(nx, nx, nc);
		WX = zeros(nx, nc);
		WXTX = zeros(nc, 1);  
		WXY = zeros(nx,ny,nc);
		
		for(cI in 1:nc) #fill (weighted by posteriors) conditional second moments. That is for each ci, Z Z^2 ZTZ and the XY  . 
		{
			#W*** seems to be wegithed moments matrices for each cluster
			#they are passed to suited functions to estimate parameters
			weights = repmat(posteriors[cI,], ny, 1); #Pr(Q=i|X)
			WYbig = Y*weights;   #estimates posterior mean E(Y!Q=i)
			WYY[,,cI] = WYbig%*%t(Y); #posterior moments of y E(y^2 | Q)
			WY[,cI] = as.matrix(rowSums(WYbig)); # E(Y|Q)
			WYTY[cI] = sum(diag(t(WYbig)%*%Y)); # seems E(YTY |Q=i)
			
			weights = repmat(posteriors[cI,], nx, 1);  # the same as above for X
			WXbig = X*weights;   
			WXX[,,cI] = WXbig%*%t(X);
			WX[,cI] = as.matrix(rowSums(WXbig));
			WXTX[cI] = sum(diag(t(WXbig)%*%X));
			WXY[,,cI] = WXbig%*%t(Y); #that is E(X) E(X^2) E(XY) in each Q cluster
		}
		#M step
		#reestimates sigmaX and muX
		revalCwrX=.mixGaussMstep(w, WX, WXX, WXTX,cov_type=cov_typeX, cov_prior=cov_priorX);
		cwrOutput$muX=revalCwrX$mu; cwrOutput$SigmaX=revalCwrX$Sigma;
		rm(revalCwrX)
		for(i in 1:nc)  #verifies positive definitive sigmaX
		{
			test=.ispossemdef(cwrOutput$SigmaX[,,i])
			if(test==FALSE) stop("Error: SigmaX matrix not semidefinite positive")
		}
		if(clamp_weights==TRUE) W = cwrStart$weightsY else W=NULL;
		#reestimates Y
		clgObj=.clgMstep(w, WY, WYY, WYTY, WX, WXX, WXY, #<80>stimates parameter of conditional linear gaussian
				cov_type= cov_typeY, clamped_weights=W,cov_prior=cov_priorY);
		cwrOutput$muY=clgObj$mu; cwrOutput$SigmaY=clgObj$Sigma; cwrOutput$weightsY=clgObj$B;
		cwrOutput$priorC = w/sum(w); #reestimates cluster membership probabilities
		
		for(i in 1:nc) #check goodness of estimates (semidef. positive)
		{
			if(!.ispossemdef(cwrOutput$SigmaY[,,i])) stop("Matrix not semidefinite positive")
		}
		
		if(verbose) cat("iteration ",  num_iter, " logLik ", logLik,"\n")
		#end M step 
		num_iter =  num_iter + 1; #increase iterations
		converged = .em_converged(logLik, previous_logLik, thresh); #checks convergence
		previous_logLik = logLik; #save convergence
		cwrOutput$posteriors=posteriors;  #print out posteriors
		
		#saves group membership
		group=apply(posteriors, 2, .maxIndex) #assigns group
		group=t(as.matrix(group)) #save group membership 
		
		cwrOutput$group=group;     #saves group membership
		cwrOutput$logLik=logLik;  #saves logLikelihood
	}
#calculates the number of parameters 
	nPar=numel(cwrOutput$muX)+numel(cwrOutput$SigmaX)+numel(cwrOutput$muY)+numel(cwrOutput$SigmaY)+numel(cwrOutput$priorC)+numel(cwrOutput$weightsY);
#AIC & BIC
	cwrOutput$nPar=nPar; #saves n of parameters
	cwrOutput$aic=-2*logLik+2*nPar; #saves aic
	cwrOutput$bic=-2*logLik+nPar*log(N); #saves bic
	cwrOutput$N=N;
#finishing
	cwrOutput$call=match.call();      #saves the call
	cwrOutput$timeTotal=(proc.time()-pmt)[3]; #saves computation time
	
	cwrOutput$X=t(X);
	cwrOutput$Y=t(Y);
	class(cwrOutput)="cwrObj"
	return(cwrOutput)
}





stepCwr<-function(X, Y, nc, prop=0.1, nIter=10, changeTrainingSet=FALSE)
{
	#force x, y to be vectors
	X=as.matrix(X);   #converts the data into matrix
	Y=as.matrix(Y);   
	if(size(X,1)!=size(Y,1)) stop ("Number of observations differs by X and Y")
	#sample the data
	howMany=round(size(X,1)*prop)        #sample a % of data
	index=sample(1:size(X,1), howMany)
	X4test=X[index,];     #sample the data
	Y4test=Y[index,];
	#initial settngs   of the EM restarter
	llCompare=-Inf;
	cwrToGive=NULL;    #this is the initial cwrObject fitted
	for(i in 1:nIter)
	{
		if(changeTrainingSet==TRUE) #if training set are needed to change
		{                         #changes each iteration tre training set
			index=sample(1:size(X,1), howMany) #to obtain best parameters.
			X4test=X[index,];    #reselect the parameters
			Y4test=Y[index,];
		}
		cwr4Test=NULL;
		cwr4Test=try(cwrEm(X4test, Y4test, nc, verbose=FALSE));      #try needed to avoid stop due to numerical errors
		if(class(cwr4Test)=="cwrObj") {
			if(cwr4Test$logLik>llCompare) #if better llikelihood
			{ #saves the output
				cwrToGive=cwr4Test; 
				llCompare=cwr4Test$logLik;
			}}
		cat("iteration ",i, " maxLL ", llCompare, "\n")
	}
	
	finalCwrObj=cwrEm(X,Y, nc,    #sends the cwrobject with inital parameters to the final fitting
			create_init_params=FALSE, cwrStart=cwrToGive, verbose=TRUE)
	return(finalCwrObj) #return
}


bestPermutation<-function(origClass, inizOutput) #orig: valori veri, inizOutput: classificazione
{
	#require(vegan)
	out=list() 
	inizOutput=as.numeric(inizOutput) #conversione necessaria
	.Sostituisci=function(cosa, conCosa) #cosa è in vettore di numeri indice di gruppi
	{                                    #conCosa è una permutazione di questi numeri indice
		originali=sort(unique(cosa))
		if(length(originali)!=length(conCosa)) stop("Errore! Numero di livelli incongruente")
		.swap=function(x, subVector) {return(subVector[x])} #uso il vettoriale x velocita   
		out=sapply(cosa, .swap, conCosa) #x=cosa , conCosa = subVector
		return(out)
	}
	
	.accuracy=function(classVal, trueVal) #restituisce accuratezzaClassificazione
	{  tabella=table(classVal, trueVal) #contingency table
		out=sum(diag(tabella)/sum(rowSums(tabella))) #diagonal proportion
		return(out)
	}
	
	nGroups=max(unique(inizOutput)) #i gruppi sono assunti indicagi con label numerii da 1,2,.,k
	
	if(nGroups==1) {
		warning("Warning! No more than one group classified"); out$permutazione=1; 
		out$group=inizOutput;
		return(out) }
	
	permutazione=1:nGroups #questa e' la classificazione originaria 1,2,.,k (nessuno scambio)
	accuratezza=.accuracy(origClass, inizOutput)
	
	
	
	matricePermutazioni=allPerms(nGroups)
	
	if(nGroups==2) matricePermutazioni=t(as.matrix(matricePermutazioni)) #problemaConversione quando nGroups==2
	
	
	for(i in 1:nrow(matricePermutazioni))
	{
		permTest=matricePermutazioni[i,] #crea la nuova permutazione delle classificazioni
		newOutput=.Sostituisci(inizOutput, permTest) #scambia la classificazione in coerenza con la nuova
		accTest=.accuracy(origClass, newOutput) #misura la classificazione
		if(accTest>accuratezza) #se e' migliore la scambia
		{
			accuratezza=accTest
			permutazione=permTest
		}
	}
	
	#sceglie la classificazione migliore (puo' essere quella di partenza
	miglioreClassificazione=.Sostituisci(inizOutput, permutazione)
	
	out$permutation=permutazione #usato piu' in la
	out$groups=miglioreClassificazione
	class(out)="bestPerm"
	return(out)
}