# TODO: Add comment
# 
# Author: Giorgio Alfredo Spedicato
###############################################################################


#makes trace of a matrix
.traceMatr<-function(A)
{
	return(sum(diag(A)))
}

.swapperCwr<-function(cwrObj, bestPermObj)
{
	origCwr=cwrObj
	nGroups=nrow(origCwr$priorC) #takes numer of groups
	#swaps mean and variances
	cwrObj$muX[,1:nGroups]=origCwr$muX[,bestPermObj$permutation]
	cwrObj$muY[,1:nGroups]=origCwr$muY[,bestPermObj$permutation]
	cwrObj$SigmaX[,,1:nGroups]=origCwr$SigmaX[,,bestPermObj$permutation]
	cwrObj$SigmaY[,,1:nGroups]=origCwr$SigmaY[,,bestPermObj$permutation]
	#swaps weights and priors
	cwrObj$weightsY[,,1:nGroups]=cwrObj$weightsY[,,bestPermObj$permutation]
	cwrObj$priorC[1:nGroups,]=cwrObj$priorC[bestPermObj$permutation,]
	#swaps posteriors
	cwrObj$posteriors[1:nGroups,]=cwrObj$posteriors[bestPermObj$permutation,]
	#corrects groups
	cwrObj$group=as.matrix(bestPermObj$groups)
	return(cwrObj)
}

# .sqdist      Squared Euclidean or Mahalanobis distance.

.sqdist<-function(p, qV, A=NULL)
{

# .sqdist(p,q)   returns m(i,j) = (p(:,i) - q(:,j))'*(p(:,i) - q(:,j)).
# .sqdist(p,q,A) returns m(i,j) = (p(:,i) - q(:,j))'*A*(p(:,i) - q(:,j)).
	
#  From Tom Minka's lightspeed toolbox
	
	d=size(p)[1]; pn=size(p)[2];
	d=size(qV)[1]; qn=size(qV)[2];
	
	if(is.null(A))
	{
		pmag = colSums(p*p);
		qmag = colSums(qV*qV);
		m = repmat(qmag, pn, 1) + repmat(pmag, 1, qn) - 2*t(p)%*%qV; #necessaria una modifica
		#m = repmat(qmag,1,pn) + repmat(t(pmag), 1,qn) - 2*t(p)%*%qV;
		#m = ones(pn,1)*qmag + pmag'*ones(1,qn) - 2*p'*q;
	} else {
		if(isempty(A)|isempty(p)) stop("Error!.sqdist: empty matrices");
		Ap = A%*%p;
		Aq = A%*%qV;
		pmag = colSums(p*Ap);
		qmag = colSums(qV*Aq);
		m = repmat(qmag, pn, 1) + repmat(pmag, 1, qn) - 2*t(p)%*%Aq;
	}
	return(m)
}

#funzione mix>GaussMstep
.mixGaussMstep<-function(w, Y, YY, YTY, cov_type="full", cov_prior=NULL, clamped_mean=NULL, clamped_cov=NULL, tied_cov=FALSE)
{
	moments=list()        #moments list
	Ysz=size(Y)[1]; Q = size(Y)[2];   #Ysz: number of elements in Y; Q=number of cluster
	N = sum(w);                       #number of observation
	
	if(is.null(cov_prior)==TRUE) #creates an ipotetical (spherical) covariance matrix
	{
		cov_prior = repmat(0.01*eye(Ysz,Ysz),  1, 1, Q);
		if(Q==1) #correction if only one cluster
		{
			tempMatr=0.01*eye(Ysz,Ysz)
			cov_prior=array(0, c(Ysz, Ysz, Q))
			cov_prior[,,1]=tempMatr 
		}
	}
	YY = matlab::reshape(YY, c(Ysz, Ysz, Q)); #FIX: i dont' know what it does
	w = w + as.matrix(as.numeric(w==0));
	
	if (!is.null(clamped_mean))  mu = clamped_mean  else {
		mu = zeros(Ysz, Q);   
		for(i in 1:Q) mu[,i] = Y[,i] / w[i]; #estimates mean: mean has already been given in Y[,i] but must be scaled by w[i]
	}
	
	if (!is.null(clamped_cov))    #if cov given.
	{
		Sigma = clamped_cov;
		moments$mu=mu;
		moments$Sigma=Sigma;
		return(moments)
	}
	#if clamped cov==NULL
	if (tied_cov==FALSE)  #if covariance not tied
	{
		Sigma = zeros(Ysz,Ysz,Q); #initialize covariance 
		for(i in 1:Q)
		{
			if(cov_type == "spherical") #equal variances in every dimension of Y (but they can differ by cluster: not tied)
			{
				s2 = (1/Ysz)*((YTY[i]/w[i]) - t(mu[,i])%*%mu[,i] ) #calculate the variance as Y^2 (rescaled by w[i])-muY^2
				Sigma[,,i] = as.numeric(s2)*eye(Ysz); #saves variance
			} else {
				SS = YY[,,i]/w[i]  - mu[,i]%*%t(mu[,i]); #generic case
				if(cov_type=="diagonal") SS = diag(diag(SS)) #if covariance is ellipsoidal
				Sigma[,,i] = SS;   #saves for each cluster
			}
		}
	} else   #tied
	{
		if(cov_type=="spherical") #
		{
			s2 = (1/(N*Ysz))*(sum(YTY,2) + sum(diag(t(mu)*mu)*w));   #this is equation 19
			Sigma = as.numeric(s2)*eye(Ysz); #s2*identityMatrix
		} else {    #if it is diagonal or full
			SS = zeros(Ysz, Ysz);
			for(i in 1:Q)
				SS = SS + YY[,,i]/N - mu[,i]%*%t(mu[,i]);   #this is equation 15
		}
		if(cov_type=="diagonal") Sigma = diag(diag(SS)) else Sigma = SS;    #if it is diagonal takes drops out the off - diagonal elements
	}
	if(tied_cov==TRUE)  Sigma =  repmat(Sigma, c(1, 1, Q));   #tied means idential sigma for each cluster / group
	Sigma = Sigma + cov_prior;
	moments$mu=mu;
	moments$Sigma=Sigma;
	return(moments)
}

.mixgauss_prob<-function(dataMatr, mu, Sigma, mixmat=NULL, unit_norm=NULL)
{
	if((min(size(mu))==1) & (size(mu,2)==1)) #verifica se è un vettore e un gruppo
	{
		d=length(mu);        
		Q=1; M=1;
	} else {
		if(ndims(mu)==2)       #execute some controls
		{
			d=size(mu)[1];
			Q=size(mu)[2];
			M=1;
		} else           
		{
			d=size(mu)[1];
			Q=size(mu)[2];
			M=size(mu)[3];
		}
		
	}
	d=size(dataMatr)[1]; T=size(dataMatr)[2]; #T=nobs , d=dimension of X
	
	#if nargin < 4, mixmat = ones(Q,1); end #io ipotizzo sempre che gli vengano passati 3 argomenti
	#if nargin < 5, unit_norm = 0; end
	mixmat=ones(Q,1); unit_norm=FALSE;
	if(sum(size(Sigma))==3) #if isscalar(Sigma) #correggo adattando per la struttura di sigma
	{
		mu=matlab::reshape(mu, d, Q*M)
		if(unit_norm==TRUE)
		{
			cat("unit norm","\n")
			Di = 2 - 2*(t(mu)%*%dataMatr);
			
			D2 = t(.sqdist(dataMatr, mu)); 
			#assert(approxeq(D,D2))
		} else {
			Di=t(.sqdist(dataMatr, mu));
		}
		rm(mu, dataMatr)
		logB2 = -(d/2)*log(2*pi*Sigma) - (1/(2*Sigma))*Di; # det(sigma*I) = sigma^d
		B2 = matlab::reshape(exp(logB2), c(Q,M,T));
		rm(logB2) # ATB: clear big old data
	} else {
		if(ndims(Sigma)==2) #non capisco cosa serva ma nn importa
		{
			mu=matlab::reshape(mu, c(d, Q*M));
			Di=t(.sqdist(dataMatr, mu, solve(Sigma)));
			logB2 = -(d/2)*log(2*pi) - 0.5*log(det(Sigma)) - 0.5*Di;
			#denom = sqrt(det(2*pi*Sigma));
			#numer = exp(-0.5 * D);
			#B2 = numer/denom;
			B2 = matlab::reshape(exp(logB2), c(Q, M, T));
		} else {
			if(ndims(Sigma)==3)            #questo ciclo e' quello che verra' eseguito piu' spesso
			{
				B2=zeros(Q,M,T);  
				for(j in 1:Q)
				{
					if(.isposdef(Sigma[,,j]))
					{
#                      Di = t(.sqdist(dataMatr, permute(mu[,j,], c(1,3,2)), inv(Sigma[,,j]))); #cosi' non riesco a farla girare: grave problema
#pero' senza permute va lo stesso
#dovremo vedere caso per caso
#
						Di = t(.sqdist(dataMatr, as.matrix(mu[,j]),  solve(Sigma[,,j]))); #con as matrix lo prende come colonna
						logB2 = -(d/2)*log(2*pi) - 0.5*log(det(as.matrix(Sigma[,,j]))) - 0.5*Di; #ma andra; verificare come matlab ed R trattano sottrazione di array
						B2[j,,] = exp(logB2);
					}  else stop("Error: Sigma not positive def")
					
				}
				
			} else    #caso con piu' dimensioni di tre? non lo immagino
			{
				B2=zeros(Q,M,T)       #queste sono le probabilita' condizionate
				for(j in 1:Q)
				{
					for(k in 1:M)
					{
						B2[j,k,]=.gaussian_prob(dataMatr, mu[,j,k], Sigma[,,j,k]) #€stimates gaussian probabilities
					}
				}
			}
		}
	}
	
	B = zeros(Q,T);
	if(Q < T)
	{
#  for (qI in 1:Q)  B[qI,] = mixmat[qI,]%*%permute(B2(qI,,), c(2, 3, 1)); #non implementato
		for (qI in 1:Q)   B[qI,] = mixmat[qI,]%*%t(as.matrix(B2[qI,,]))
		#B(q,:) = mixmat(q,:) * squeeze(B2(q,:,:)); % squeeze chnages order if M=1
		
	} else {
		for (tI in 1:T)
		{B[,tI] = rowSums(mixmat*B2[,,tI])}
	}
	mixgaussObj=list()
	mixgaussObj$B=B
	mixgaussObj$B2=B2
	return(mixgaussObj)
}

#programma mixgaus init
.mixgauss_init<-function(M, dataM, cov_type="full", method="kmeans")  #nc, X, cov_typeX,.
{
	d=size(dataM)[1]   #number of dimension
	T=size(dataM)[2]   #number of obs
#data = reshape(data, d, T); #boh?
	Sigma=array(0, c(d,d, M)) #dxd is the var covariance matrix for each of M dimensions
	
	if(method=="rnd")
	{
		Ci=cov(t(dataM))   #create a covariance matrix from given data
		#Sigma = repmat(diag(diag(C))*0.5, [1 1 M]);
		SigmaIntermedia=diag(diag(Ci))*0.5 #takes only diagonal (seems)
		for(i in 1:M) Sigma[,,i]=SigmaIntermedia  #matrice diagonale in ogni sottogruppi
		mu=dataM[,sample(1:T, M)] #prende due punti qualsiasi  per la media
		pesi=rep(1,M)
		pesi=pesi/sum(pesi) #equal group weights
		weights=as.matrix(pesi) #creates a vector
	}
	if(method=="kmeans")
	{
		#require(cluster)
		#c'è un problema: implemento solo le le matrici
		#con varianza covarianza completa
		kmeans4dataM=kmeans(t(dataM),M) #esegue kmenas   per 2 gruppi
		dataMRielaborato=as.data.frame(t(dataM)); #trasforma temporaneamente in df
		cluId=as.numeric(kmeans4dataM$cluster); #prende il clusterId
		mu=matrix(0,d,M)   #mean matrix
		pesi=numeric(0)
		for(i in 1:M)
		{
			indicizzatore=which(cluId==i) 
			dataTemp=as.matrix(dataMRielaborato[indicizzatore,]) #selects only i-th cluster elements
			meanI=as.matrix(apply(dataTemp,2,mean))  
			mu[,i]=meanI #calculates means and savels them
			SigmaIntermedia=if(min(size(dataTemp))==1) var(dataTemp) else cov(dataTemp) 
			SigmaIntermedia[is.na(SigmaIntermedia)]=1 #calculates var / covar matrices
			Sigma[,,i]=SigmaIntermedia
			lunClu=dim(dataTemp)[1] 
			pesi=c(pesi, lunClu)    #conunts number of elements
		}
		pesi=pesi/sum(pesi) #saves weights
		weights=as.matrix(pesi)
	}
	initOutput=list() #returns initial ouputpu
	initOutput$mu=mu
	initOutput$weights=weights
	initOutput$Sigma=Sigma
	return(initOutput)
}

#funzione is pos def
.isposdef<-function(Matr)
{
	out=FALSE;
	if(all(eigen(Matr,only.values=T)$values>0)) out=TRUE
	return(out)
}

#funzione is pos semi def
.ispossemdef<-function(Matr)
{
	out=FALSE;
	if(all(eigen(Matr,only.values=T)$values>=0)) out=TRUE
	return(out)
}

#funzione x la probabuilità gaussiana
.gaussian_prob<-function(x, m, C, use_log=FALSE)
{
	if(max(dim(m))==1) x=t(c(x))
	d=size(x)[1];
	N=size(x)[2];
	m = t(c(m));
	M = m*ones(1,N);
	denom = (2*pi)^(d/2)*sqrt(abs(det(C)));
	mahal = apply((t(x-M)%*%solve(C)*t(x-M)),1, sum);
	if(any(mahal<0))   warning("mahal < 0 => C is not psd")
	if (use_log==TRUE) p = -0.5*mahal - log(denom) else p = exp(-0.5*mahal) / (denom+0.00000000000001)
	return(p)
}

#funzione x verificare la confergenza
.em_converged<-function(logLik, previous_logLik, thresh, check_increased=TRUE)
{
	
	if(!is.finite(previous_logLik)) return(FALSE) #alla prima iterazione è -infiito
	converged = FALSE;
	decrease = FALSE; #parametri iniziali
	
	if(check_increased) #verifies if logLikelihood has decreased
	{
		if(logLik - previous_logLik < 0)
		{
			cat("likelihood decreased from ", previous_logLik," to ", logLik,"\n");
			decrease = TRUE;
			converged = FALSE;
			return(converged);
		}
	}
	delta_logLik = abs(logLik - previous_logLik); #has converged if deltalik/pastlik < thrshold
	avg_logLik = (abs(logLik) + abs(previous_logLik))/2;
	ratio=(delta_logLik / avg_logLik)
#cat("delta ",delta_logLik, " avg ", avg_logLik, " ratio ", ratio, "\n")
	if(ratio<thresh) converged = TRUE;
	return(converged)
}

#funzione cwr prob
.cwr_prob<-function(cwrObj, X, Y)   #seems to work!
{
	probList=list()   #initial creation of the list
	nx=size(X,1);    #how many x vars
	N=size(X,2);     #how many obs
	nc=length(cwrObj$priorC) #how many clusters
	
	if (nc==1)                            #should be checked
	{
		predCwr=.cwr_predict(cwrObj,X);
		mu=predCwr$mu;
		Sigma=predCwr$Sigma;
		probList$likY=.gaussian_prob(Y, mu, Sigma);
		probList$likXandY=likY;     #X do no exist
		probList$likYgivenX = likY;
		probList$post = rep(1,N); #equal probability
		
	} else {
		#this part seems to estimate probabilities 
		#check carefully 
		likY = .clg_prob(X, Y, cwrObj$muY, cwrObj$SigmaY, cwrObj$weightsY); #calculates P(Y|Q=i)
		likXObj = .mixgauss_prob(X, cwrObj$muX, cwrObj$SigmaX) #check what it does?
		likX=likXObj$B2  #posterior probabilities
		likX = drop(likX); #calculates P(X|Q=i) in matrix form
		prior = repmat(cwrObj$priorC, 1, N); #verificare! possibile incongruenza
		
		postProb = likX*likY*prior; #why posterior probabilies seems to be indipendent?
		likXandY = t(as.array(colSums(postProb))); 
		postProb = postProb/repmat(t(likXandY), nc, 1); #these are posterior probabilities
		
		likX = colSums(likX*prior); #pr(X) non contizionata
		likYgivenX = likXandY/likX;   #conditional probability
		probList$likX=likX; #saves and returns probability
		probList$likYgivenX=likYgivenX ;
		probList$post=postProb;
		probList$likY=likY;
		probList$likXandY=likXandY;
	}
	return(probList);
}

#funzione cwr predict
.cwr_predict<-function(cwrObj, X) #non ho implementato l'argomento mask
{
	cwrPred=list()
	dimX=size(X); nx=dimX[1]; T=dimX[2];
	dimWeightsY=size(cwrObj$weightsY);ny=dimWeightsY[1];nx=dimWeightsY[2];nc=dimWeightsY[3]
	mu = zeros(ny, T);
	igma = zeros(ny, ny, T);
	
#pongo comp_mask=false
	comp_mask=FALSE
	
	if (nc==1)
	{
		if(is.null(cwrObj$weightsY)) #copio pari pari da matlab
		{
			mu=repmat(cwrObj$muY, 1, T)
			Sigma = repmat(cwrObj$SigmaY, 1, 1, T); # e' una prova ma sembra averlo preso
		} else {
			mu = repmat(cwrObj$muY, 1, T) + cwrObj$weightsY%*%X; #anche questa e' una prova
			Sigma = repmat(cwrObj$SigmaY, 1, 1, T);
		}
		
		cwrPred$mu=mu; #qua ci sarebbe un po' di mask non implementato
		cwrPred$Sigma=Sigma;
		cwrPred$weights=weights;
		return(cwrPred)
	} else { #se i gruppi sono + di uno
		likObj=.mixgauss_prob(X, cwrObj$muX, cwrObj$SigmaX); #implementare .mixgauss_prob
		#weights = normalize(repmat(cwr.priorC, 1, T) .* likX, 1);
		#devono sommare per colonna ad uno
		weights = repmat(cwrObj$priorC, 1, T)*likObj$likX;
		sColWeights=apply(weights , 2, sum)
		weights=weights / sColWeights;
		
		for (t in 1:T)
		{
			mut=zeros(ny,nc);
			for (c in 1:nc)
			{
				mut[,c]=cwrObj$muY[,c]+cwrObj$weightsY[,,c]%*%X[,t];
				#non e' stato implementato comp_mask
			}
			collapseObj=.collapse_mog(mut, cwrObj$SigmaY, cwrObj$weights[,t]) #implementare e verificare
			mu[,t]=collapseObj$new_mu
			Sigma[,,t]=collapseObj$Sigma
		}
		
		cwrPred$mu=mu; #qua ci sarebbe un po' di mask non implementato
		cwrPred$Sigma=Sigma;
		cwrPred$weights=weights;
		return(cwrPred)
		
	}
	
}


.collapse_mog<-function(mu, Sigma, coefs)
{
	collapseObj=list() #inizializzo la lista
	new_mu=apply(mu%*%diag(coefs), sum, 2)
	n = length(new_mu);
	new_Sigma = zeros(n,n);
	new_Sigma2 = zeros(n,n);
	
	for (j in 1:length(coefs))
	{
		m=mu[,j]-new_mu;
		new_Sigma = new_Sigma + coefs[j]%*%(Sigma[,,j] + m%*%t(m));
		new_Sigma2 = new_Sigma2 + coefs[j]%*%(Sigma[,,j] + mu[,j]*t(mu[,j]));
	}
	
	collapseObj$new_mu=new_mu;
	collapseObj$new_Sigma=new_Sigma;
	collapseObj$new_Sigma2=new_Sigma2;
}

#internal function to find which group (column) has max posterior probability
.maxIndex=function(colVec)  {
	return(which(colVec==max(colVec)))
} 


.clgMstep<-function(w, Y, YY, YTY, X, XX, XY,cov_type= "full", tied_cov=FALSE,  clamped_cov=NULL, clamped_mean=NULL,clamped_weights=NULL, cov_prior=NULL, xs=NULL, ys=NULL, post=NULL)
{
	clgObj=list()
	Ysz=size(Y)[1]; Q = size(Y)[2];   #Ysz number of Y dim (variable) Q number of clusters
	
	if(missing(X)) #lo pongo missing: no regression
	{
		#%B = [];
		#B2 = zeros(Ysz, 1, Q);
		#B=zeros(Ysz, 1, Q); #patch
		#for (i in 1:Q) B[,,i] = zeros(Ysz,1) #verificare
		momentsObj = .mixGaussMstep(w, Y, YY, YTY, cov_type="full", cov_prior=cov_prior, clamped_mean=clamped_mean, clamped_cov=clamped_cov, tied_cov=tied_cov)   #ristima i momenti di (solo Y)
		clgObj$mu=momentsObj$mu; clgObj$Sigma=momentsObj$Sigma; clgObj$B=NULL;
		return(clgObj) #restituisce il vettore
	} #non sara' testato per ora
	
	N = sum(w);    #number ob observations
	if(is.null(cov_prior)==TRUE)
	{
		cov_prior = 0.01*repmat(eye(Ysz,Ysz), c(1, 1, Q));
		if(Q==1)   #force array in R when Q=1: correction needed if Q=1
		{
			tempMatr=0.01*eye(Ysz,Ysz)
			cov_prior=array(0, c(Ysz, Ysz, Q))
			cov_prior[,,1]=tempMatr 
		}  
	}
	w = w +as.matrix(as.numeric(w==0)); # set any zero weights to one before dividing: numerical (?) correction
	Xsz = size(X,1);  #number of X dimension
	ZZ = zeros(Xsz+1, Xsz+1, Q); #as written in Murphy article a one is added to   ZZ = E((X, 1)(X 1)')(
	ZY = zeros(Xsz+1, Ysz, Q);  #ZY = (X 1)Y'
	
	if(Q==1)  #r creates matrices in place of array when third dimension is 1
	{
		#converts XX in array
		tempDim=as.numeric(size(XX));
		copyTemp=array(0, c(tempDim, 1));
		copyTemp[,,1]=XX;
		XX=copyTemp;
		#converts XY in array
		tempDim=as.numeric(size(XY));
		copyTemp=array(0, c(tempDim, 1));
		copyTemp[,,1]=XY;
		XY=copyTemp;
		#converts YY in array
		tempDim=as.numeric(size(YY));
		copyTemp=array(0, c(tempDim, 1));
		copyTemp[,,1]=YY;
		YY=copyTemp;
		#converts X in array
		#tempDim=as.numeric(size(X));
		#copyTemp=array(0, c(tempDim, 1));
		#copyTemp[,1]=X;
		#X=copyTemp;
		#converts Y in array
		#tempDim=as.numeric(size(Y));
		#copyTemp=array(0, c(tempDim, 1));
		#copyTemp[,1]=Y;
		#Y=copyTemp;
	}
	
#sono tutte equazioni articolo murphy clg
#does for eacgh cluster 
	for(i in 1:Q)
	{
		ZZ[,,i] = rbind(cbind(XX[,,i],  as.matrix(X[,i])),cbind(t(X[,i]),w[i])); #all present in eq 9  (up)
		ZY[,,i] = rbind(as.matrix(XY[,,i]),t(as.matrix(Y[,i])));  #see eq 9  (down)
		#correzioni di conversione fixed: old rbind(XY[,,i],Y[,i]); correzione 2 as.matrix(Y[,i]) senza trasposto
	}
	
#vari casi
	
	if ((!is.null(clamped_weights)) & (!is.null(clamped_mean))) 
	{
		B = clamped_weights; #if are passed
		mu = clamped_mean;
	}
	
	if (!is.null(clamped_weights) & is.null(clamped_mean))
	{
		B = clamped_weights; #if only clamped_weights are passed
		# this is quation n5 5
		mu = zeros(Ysz, Q); #initialize mena
		for(i in 1:Q) mu[,i] = (Y[,i] - B[,,i]%*%X[,i]) / w[i]; #formula 7
	}
	
	if(is.null(clamped_weights) & !is.null(clamped_mean)) #testare
	{
		mu = clamped_mean;
		# this is qeuation n  3
		B = zeros(Ysz, Xsz, Q);
		for(i in 1:Q)
		{
			tmp = t(XY[,,i]) - mu[,i]*t(X[,i]);   #see equation n 3
			#B(:,:,i) = tmp * inv(XX(:,:,i));
			B[,,i] = t(solve(XX[,,i] , t(tmp)));   #see equation n 3
		}
	}
	
	if(is.null(clamped_weights) & is.null(clamped_mean))
	{
		mu = zeros(Ysz, Q); #iniktialize as zero
		B = zeros(Ysz, Xsz, Q);
		# Nothing is clamped, so we must estimate B and mu jointly
		for(i in 1:Q)
		{
			#l'ho lasciato non commentato
			# eqn 9
			#if rcond(ZZ(:,:,i)) < 1e-10
			#  sprintf('clg_Mstep warning: ZZ(:,:,%d) is ill-conditioned', i);
			# probably because there are too few cases for a high-dimensional input
			# ZZ(:,:,i) = ZZ(:,:,i) + 1e-5*eye(Xsz+1);
			#end
			#%A = ZY(:,:,i)' * inv(ZZ(:,:,i));
			A = t(solve(ZZ[,,i] , as.matrix(ZY[,,i]))); # FIX    A = t(solve(ZZ[,,i] , ZY[,,i]));   #this is equation 9
			B[,,i] = A[, 1:Xsz];     #B, weights are the first Xsz columns and mu (intercept) are the last column
			mu[,i] = A[, Xsz+1];
		}
		
	}
	
	if(!is.null(clamped_cov)) #se la cov. viene data allora
	{                         #restituisce quanto elaborato (B, mu, Sigma)
		Sigma = clamped_cov;
		clgObj$Sigma=Sigma;
		clgObj$mu=mu;
		clgObj$B=B;
		return(clgObj)
	}
	
	if(cov_type=="spherical")
	{
		if(tied_cov==FALSE)
		{
			Sigma=zeros(Ysz, Ysz, Q);
			for(i in 1:Q)
			{
				# this is eqn 16
				A = cbind(t(as.matrix(B[,,i])), as.matrix(mu[,i])); 
				#A = cbind(B[,,i], mu[,i]); #correggo ancora!: A = cbind(t(as.matrix(B[,,i])), t(as.matrix(mu[,i])))
				#s = trace(YTY(i) + A'*A*ZZ(:,:,i) - 2*A*ZY(:,:,i)) / (Ysz*w(i)); % wrong!
				s = (YTY[i] + .traceMatr(t(A)%*%A%*%ZZ[,,i]) - .traceMatr(2*A%*%ZY[,,i])) / (Ysz*w[i]);
				Sigma[,,i] = s*eye(Ysz,Ysz);    #does this work for alla cluster
				#ho lasciato tale e quale il commentro trovato in Murphy
				#%%%%%%%%%%%%%%%%%%% debug
				#if ~isempty(xs)
				#[nx T] = size(xs);
				#zs = [xs; ones(1,T)];
				#yty = 0;
				#zAAz = 0;
				#yAz = 0;
				#	for t=1:T
				#yty = yty + ys(:,t)'*ys(:,t) * post(i,t);
				#zAAz = zAAz + zs(:,t)'*A'*A*zs(:,t)*post(i,t);
				#yAz = yAz + ys(:,t)'*A*zs(:,t)*post(i,t);
				#end
				#assert(approxeq(yty, YTY(i)))
				#assert(approxeq(zAAz, trace(A'*A*ZZ(:,:,i))))
				#assert(approxeq(yAz, trace(A*ZY(:,:,i))))
				#s2 = (yty + zAAz - 2*yAz) / (Ysz*w(i));
				#	assert(approxeq(s,s2))
				#end
				#%%%%%%%%%%%%%%% end debug
				
			}
		}  else {       #no tied
			S = 0;
			for(i in 1:Q)
			{
				# this is equation  eqn 18
				A = cbind(t(as.matrix(B[,,i])), as.matrix(mu[,i])); #A = cbind(B[,,i], mu[,i]); #correggo ancora!: A = cbind(t(as.matrix(B[,,i])), t(as.matrix(mu[,i])))
				S = S + .traceMatr(YTY[i] + t(A)%*%A%*%ZZ[,,i] - 2*A%*%ZY[,,i]); #almeno nell'esempio che ho implementato
			}
			#the estimated sigma is replied for all clusters
			Sigma = repmat(S / (N*Ysz), c(1, 1,Q));
		}
	} else {        #full diagonal
		if(tied_cov==FALSE)
		{
			Sigma = zeros(Ysz, Ysz, Q);
			for(i in 1:Q)
			{
				A = cbind(t(as.matrix(B[,,i])), as.matrix(mu[,i])); #A = cbind(B[,,i], mu[,i]); #correggo ancora!: A = cbind(t(as.matrix(B[,,i])), t(as.matrix(mu[,i])))
				# eqn 10
				SS = (YY[,,i] - t(ZY[,,i])%*%t(A) - A%*%ZY[,,i] + A%*%ZZ[,,i]%*%t(A)) / w[i];
				if(cov_type=="diagonal") Sigma[,,i] = diag(diag(SS))
				else Sigma[,,i] = SS;
			}
		} else { #tied
			SS = zeros(Ysz, Ysz);
			for(i in 1:Q)
			{
				A = cbind(t(as.matrix(B[,,i])), as.matrix(mu[,i])); #A = cbind(B[,,i], mu[,i]); #correggo ancora!: A = cbind(t(as.matrix(B[,,i])), t(as.matrix(mu[,i])))
				# equazione 13
				SS = SS + (YY[,,i] - t(ZY[,,i])%*%t(A) - A%*%ZY[,,i] + A%*%ZZ[,,i]%*%t(A)); #cal
			}
			SS = SS / N;  #this is "mean variance matrix"
			if(cov_type=="diagonal") Sigma = diag(diag(SS)) else Sigma = SS;
			Sigma = repmat(Sigma, c(1, 1, Q));# replies for all cluster
		}
		Sigma=Sigma+cov_prior;  #adds 0.01 diagonal
	}
	clgObj$mu=mu;
	clgObj$Sigma=Sigma;
	clgObj$B=B;
	return(clgObj)
} #fine funzione

.clg_prob<-function(X, Y, mu, Sigma, W) #it works (pare)
{
	d=size(Y)[1]; T=size(Y)[2]; #d=number of dimension, T number of samples
	d=size(mu)[1]; nc=size(mu)[2]; #overrides above
	p=zeros(nc, T) #initial sets 
	for (k in 1:nc)
	{
		denom=(2*pi)^(d/2)*sqrt(abs(det(as.matrix(Sigma[,,k]))));  #this is the denominator of a multivariate normal
		M=repmat(mu[,k], 1, T)+W[,,k]%*%X; #replies the local mean vector for all T sample elements
		sMenoUno=solve(Sigma[,,k]); #inverts Sigma
		mahal = apply((t(Y-M)%*%sMenoUno)*t(Y-M),1, sum); #evaluate args in exp(-.5*args) in the evaluation of posterior probability
		p[k,] = t(exp(-0.5*mahal) / denom);
	}
	return(p)
}
