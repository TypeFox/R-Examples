#' Primary function for hierarchical, areal analysis of distance sampling data.  This function
#' pre-processes data and calls other functions to perform the analysis, and is the only function
#' the user needs to call themselves. 
#'
#' @param Dat 	A data frame with the following columns:
#' 		(1)transect ID; 
#' 		(2)match number  currently, a maximum of 2 observers on each transect;
#' 		(3)(Observer ID);
#' 		(4)(Observation (0/1));
#' 		(5) Observed species (integer - the max integer being 'unknown' if applicable) [NOTE: modeled as factor, but need to be input as integers to account for unknown species observations]
#' 		(6-x)(Observer covariates); (things like survey conditions or observer skill; things that don't change during a transect.  Note these also need to be provided in Obs.cov)
#' 		(x+1)(Distance; if all integers, assumed to be discrete bins; if continuous, assumed standardized to (0,1) interval);
#' 		(x+2-??)(Group size and other individual covariates thought to influence detection; if group size is one of them, it's assumed to be column x+2);
#' 		Note that column names can be used to tag covariates, and that object types (e.g. numeric, factor) will be preserved in analysis
#' @param Adj   Adjacency matrix for habitat cells (diagonal matrix implies spatial independence)
#' @param Area.hab   A vector giving the area of each geographical strata (default is equal area)
#' @param Mapping  A vector giving the habitat cell id number for each transect
#' @param Area.trans	A vector giving the effective area covered by each transect as fraction of total area in the strata it is located
#' @param Observers	A (2 x number of transects) matrix giving the observers IDs that were present for each transect (the 2nd row is to contain NAs if only 1 observer was present)
#' @param Bin.length	If distances are binned, this vector gives the relative length of each distance bin (vector must sum to one) 
#' @param n.obs.cov	Number of observer covariates (e.g., seat position, visibility, etc.)
#' @param Hab.cov	A data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates
#' @param Obs.cov  A (max # of observers X # of transects X # of observer covariates) size array giving observer covariate values for each transect flown
#' @param Hab.pois.formula	A list of formulas (one for each species) giving the specific model for Poisson abundance intensity at the strata level (e.g., ~Vegetation+Latitude) for each species
#' @param Hab.bern.formula  If ZIP=TRUE, a list of formulas (one for each species) giving the specific model for the zero component for abundance intensity at the strata level (e.g., ~Vegetation+Latitude) for each species
#' @param detect If TRUE (the default), detectability is estimated; if FALSE, assumes detection probability is 1.0 (i.e. a strip transect with perfect detection).  
#' @param Det.formula  A formula giving the model for detection probability (e.g. ~Distance+Group+Visibility+Observer). Note that
#'				there are several "reserved" variable names.  "Distance", "Observer", "Species", and "Group" are reserved variable names.
#' @param Cov.prior.pdf	If individual covariates are provided, this character matrix gives the form of the prior pdfs for each covariate
#'		  current possibilities are "poisson", "pois1","poisson_ln","pois1_ln",uniform.disc","multinom","uniform.cont", or "normal".
#'		  "pois1" is 1+x where x~poisson; "poisson_ln" and "pois1_ln" are lognormal poisson models that incorporate overdispersion.  Note
#'        the dimension of this matrix are (# species X # of covariates)
#' @param Cov.prior.parms	A (s X k X n) array where s is the number of species, n is the number of individual covariates (other than distance), and
#' 		k is the maximum number of parameters considered for a single covariate (NAs can be used to fill this matrix
#'      out for covariate priors that have <k parameters).  If Cov.prior.fixed=1 for a given entry, the prior parameters supplied
#'      in each column apply to the prior pdf itself, and are treated as fixed.  If Cov.prior.fixed=0, the model will attempt
#'  	to estimate the posterior distribution of model parameters, given hyperpriors.  In this case, it is actually the hyperpriors
#'      that are being specified.  For "poisson", and "pois1", it is assumed that lambda~gamma(alpha,beta), so alpha
#' 		and beta must be supplied.  For "poisson_ln", and "pois1_ln", the model is lambda_i=exp(-sigma*Z_i+theta), so it is priors
#' 		for theta and sigma that are specified (in that order).  Theta is assumed to have a normal(mu,s^2) distribution,
#' 		and sigma is assumed to have a uniform(0,a) distribution; thus, priors are specified for these models as (mu,s, and a).
#' 		For the multinomial pdf, prior parameters of the dirichlet distribution must be specified if Cov.prior.fixed=1.
#' @param Cov.prior.fixed  An indicator matrix specifying which (if any) individual covariate distributions should be fixed during estimation
#' @param Cov.prior.n  An (# species X # indiv. covariates) matrix giving the number of parameters in each covariate pdf
#' @param pol.eff 	For continuous distance, which polynomial degrees to model (default is c(1:2); an intercept is always estimated when "Distance" is listed in "Det.formula")
#' @param ZIP  If TRUE, estimate ZIP model for abundance that includes a Bernoulli model for zeros and a Poisson + 1 model for positive values (default is FALSE)
#' @param point.ind  Estimate a correlation parameter for detection probability that's an increasing function of distance?
#' @param spat.ind	If TRUE, assumes spatial independence (no spatial random effects on abundance intensity) default is FALSE
#' @param last.ind If point independence is modeled (point.ind=TRUE), last.ind=TRUE will set observer covariance to zero for the greatest distance and maximal correlation in first bin (default is FALSE)
#' @param cor.const If TRUE, forces estimates of correlation associated with point independence to be positive if last.ind==FALSE or negative if last.ind==TRUE (default is FALSE)
#' @param fix.tau.nu  If TRUE, fixes tau.nu during estimation (the value to fix it to can be provided in "Inits")
#' @param srr  If TRUE, uses spatially retricted regression, where smoothing occurs on residuals and all spatial effects are orthogonal to the linear predictors (by default, analysis is limited to the highest 50 eigenvalues of the decomposition of the residual projection matrix to reduce computing time)
#' @param srr.tol Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation (default is 0.5)
#' @param misID If TRUE, updates species for observed animals and estimates misID parameters
#' @param misID.mat With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#' @param misID.models A formula vector providing linear model-type formulas for each positive value of misID.mat.  
#' @param misID.symm If TRUE, the constraint pi^{i|j}=pi^{j|i} is implemented; in this case, entries for pi^{j|i} are all assumed to be given a '-1' in misID.mat
#' @param grps 	If FALSE, detections are assumed to all be of individual animals
#' @param M		Matrix with species-specific rows giving maximum possible value for number of groups present in each transect (in practice just set high enough that values at M and above are never sampled during MCMC) and can be fine tuned as needed
#' @param Control	A list object including the following objects:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run; 
#'	"MH.cor": Metropolis-hastings tuning parameter for updating the correlation parameter (if point.ind==TRUE);
#'	"MH.nu": MH tuning parameters for Nu parameters (dimension = # species X # of unique strata sampled)
#'	"RJ.N": A matrix giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each species (row) and each transect (column)
#'  "iter.fix.N": Number of iterations to skip RJMCMC step at beginning of estimation (useful for cases when estimation is unstable)
#' @param Inits	An (optional) list object providing initial values for model parameters, with the following objects:
#'  "hab.pois": Initial values for habitat linear predictor parameters for poisson model;
#'  "hab.bern": If ZIP=TRUE, initial values for habitat linear predictor parameters for bernoulli zero model;
#'	"det": Initial values for detection model (includes distance, observer, env. variables, and individual covariates);
#'	"cor.par": If point.ind==TRUE, this is an initial value for the correlation parameter (which must be in (0,1));	
#'	"Nu": Gives log(lambda) for each spatial strata;
#'	"Eta.pois": If spat.ind==FALSE, spatial random effects for Poisson abundance model; one for each cell and for each species
#'  "Eta.bern": If spat.ind==FALSE & ZIP=TRUE, spatial random effects for Bernoulli abundance model; one for each cell and for each species
#'	"tau.eta.pois": If spat.ind==FALSE, precision for spatial ICAR model(s) for the Poisson component
#'  "tau.eta.bern": If spat.ind==FALSE & ZIP=TRUE, precision for spatial ICAR model(s) for the Bernoulli component
#'	"tau.nu": Precision for Nu (overdispersion relative to the Poisson distribution)
#'  One need not specify an initial value for all parameter types (if less are specified, the others are generated randomly)
#' @param adapt	If adapt==TRUE, run an additional Control$adapt number of MCMC iterations to optimize MCMC proposal distributions prior to primary MCMC
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following objects
#'	"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'  "b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'	"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
#'	"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
#'	"beta.tau": prior precision for regression coefficients (assumed Normal(0,(beta.tau*X'X)^(-1))
#' @param post.loss If TRUE, calculates observed values and posterior predictions for detection data to use with posterior predictive loss functions
#' @return returns a list with the following objecs: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis-Hastings;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' @export
#' @import Matrix
#' @keywords areal model, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn \email{paul.conn@@noaa.gov} 
#' @examples print("example analysis included in the script example_analysis.R")
hierarchical_DS<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,Observers,Bin.length,Hab.cov,Obs.cov,Hab.pois.formula,Hab.bern.formula=NULL,Det.formula,detect=TRUE,Cov.prior.pdf,Cov.prior.parms,Cov.prior.fixed,Cov.prior.n,n.obs.cov=0,pol.eff=c(1:2),ZIP=FALSE,point.ind=TRUE,spat.ind=FALSE,last.ind=FALSE,cor.const=FALSE,fix.tau.nu=FALSE,srr=TRUE,srr.tol=0.5,misID=FALSE,misID.models=NULL,misID.mat=NULL,misID.symm=TRUE,Inits=NULL,grps=FALSE,M,Control,adapt=TRUE,Prior.pars,post.loss=TRUE){
	#require(mvtnorm)
	#require(Matrix)
	#require(truncnorm)
	#require(mc2d)
	#require(MCMCpack)
	DEBUG=FALSE
	
	Adj=as.matrix(Adj)  #just in case the adjacency matrix = 1 (for 1 transect)
	S=length(Area.hab)
	n.transects=length(Area.trans)
	n.ind.cov=ncol(Dat)-(6+n.obs.cov) #number of individual covariates 
	n.species=length(unique(Dat[,"Species"]))
	if(misID==TRUE){
		n.species=nrow(misID.mat)
		n.obs.species=ncol(misID.mat)
		i.unknown=n.obs.species-n.species
	}
	
	#By no means exhaustive checking to make sure input values are internally consistent
	n.obs.max=ifelse(sum(is.na(Observers[2,]))==n.transects,1,2)
	if(n.obs.max>2)cat("\n ERROR: Current max number of observers per transect is 2\n")
	if(nrow(Control$MH.nu)!=n.species)cat("\n ERROR: Control$MH.nu does not have # rows = number of species \n")
	if(nrow(M)!=n.species)cat("\n ERROR: M does not have nrow = number of species \n")
	if(ncol(M)!=n.transects)cat("\n ERROR: M does not have ncol = number of transects \n")
  if(ncol(Control$RJ.N)!=n.transects)cat("\n ERROR: # columns of Control$RJ.N not = number of transects \n")		
  if(ZIP==TRUE & is.null(Hab.bern.formula))cat("\n ERROR: must specify Hab.bern.formula when ZIP=TRUE \n")
  if(n.species!=length(Hab.pois.formula))cat("\n ERROR: make sure length of Hab.pois.formula list equal to the number of species \n")
  if(grps==TRUE & Cov.prior.pdf[1,1]!='pois1' & Cov.prior.pdf[1,1]!='pois1_ln')cat("\n ERROR: Cov.prior.pdf for group size needs to be zero truncated (pois1 or pois1_ln)")
  #More later...	
  
  #adust M to be divisible by 2
  M[which(M%%2==1)]=M[which(M%%2==1)]+1

	if(length(unique(Dat[,3]))>1)Dat[,3]=as.factor(Dat[,3])  #convert observer to factors if not already
	if(length(unique(Dat[,5]))>1)Dat[,5]=as.factor(Dat[,5])  #convert species to factors if not already
	cur.colnames=colnames(Dat)
	cur.colnames[6+n.obs.cov]="Distance"
	if(grps==TRUE)cur.colnames[7+n.obs.cov]="Group"
	if(length(colnames(Dat))!=length(cur.colnames))cat("\n ERROR: mismatch between dimension of Dat and expected columns: check to make sure n.obs.cov, etc. correct")
	colnames(Dat)=cur.colnames
	i.binned=0
	n.bins=NULL
	if(is.factor(Dat[1,"Distance"])==1 | sum(as.numeric(Dat[,"Distance"])%%1)==0){
		i.binned=1
		n.bins=length(unique(Dat[,"Distance"]))
    Dat[,"Distance"]=as.factor(Dat[,"Distance"])
	}
  
  if(i.binned==1)if(length(Bin.length)!=n.bins)cat("\n ERROR: length of Bin.length must equal n.bins")


	Dat[Dat[,"Obs"]==0,"Species"]=Dat[Dat[,"Obs"]==1,"Species"][1]  #just set missing species to the first 'real' species
  Dat[,"Species"]=as.factor(as.character(Dat[,"Species"]))  #reestablish levels
  #convert character entries to factors
  for(i in 1:ncol(Dat))if(is.character(Dat[,i])&length(unique(Dat[,i]))>1)Dat[,i]=as.factor(Dat[,i])
  
  #for factors, determine levels, save labels, and convert to numeric
	factor.ind=sapply(Dat[1,],is.factor)
  which.factors=which(factor.ind==1)
  n.factors=sum(factor.ind)
  
  Factor.labels=vector("list",n.factors)
  for(i in 1:n.factors){
    Factor.labels[[i]]=levels(Dat[,which.factors[i]])
  }
  
	Dat.num=Dat
  #if(sum(Dat[,"Obs"]==0)>0)Dat.num[Dat[,"Obs"]==0,"Species"]=0  #if a missing obs, set species=0
	for(icol in which.factors){
		Dat.num[,icol]=as.numeric((Dat[,icol]))
	}
  Levels=vector("list",n.factors)
	for(i in 1:n.factors){
	  Levels[[i]]=sort(unique(Dat.num[,which.factors[i]]))
	}
	names(Levels)=colnames(Dat[,which.factors])
	if(misID)Levels$Species=Levels$Species[-which(Levels$Species==(n.species+1))]
  
  #update observer covariate values to reflect new factor values going from 1,2,...
  if(n.obs.cov>0){
    for(icov in 1:n.obs.cov){
      if((icov+5)%in%which.factors){
        Obs.cov[,,icov]=as.factor(as.numeric(as.factor(Obs.cov[,,icov])))
      }
    }	
  }
	
	n.Observers=apply(1-is.na(Observers),2,'sum')
	if(point.ind==TRUE & max(n.Observers)==1)cat("\n ERROR: can't have point independence when there are no transects with 2 observers \n")
	
	M=t(t(M)*n.Observers)	#actual dimension of M goes up if >1 observer 
	max.M=max(M)
	Observers.num=matrix(as.numeric(as.factor(Observers)),dim(Observers))
  
	#add an additional column for "True species" and fill
	True.sp=Dat.num[,"Species"]*Dat.num[,"Obs"]
  if(sum(True.sp==0)>0 | sum(is.na(Dat.num[,"Species"]))>0){  #if some of the species values are missing, set to value for other observer 
    True.sp.miss=which(True.sp==0 | is.na(True.sp)==1)
    for(i in 1:length(True.sp.miss)){
      if(True.sp.miss[i]%%2==1)True.sp[True.sp.miss[i]]=True.sp[True.sp.miss[i]+1]
      else True.sp[True.sp.miss[i]]=True.sp[True.sp.miss[i]-1]
    }
  }
	Obs.sp=tabulate(True.sp)[1:n.species]
	unk.ind=which(True.sp==(n.species+1))
	if(length(unk.ind)>0)True.sp[unk.ind]=sample(c(1:n.species),length(unk.ind),prob=Obs.sp,replace=TRUE)
	True.sp[which(duplicated(Dat.num[,"Match"])==TRUE)]=True.sp[which(duplicated(Dat.num[,"Match"])==TRUE)-1]
	
	if(DEBUG==TRUE)True.sp=Out$True.species #for debugging
	Dat.num=cbind(Dat.num[1:5],True.sp,Dat.num[6:ncol(Dat.num)])
	
	#Initialize data augmentation multi-d array ("Data"), parameter vectors and matrices
	Data<-array(0,dim=c(n.species,n.transects,max.M,5+n.obs.cov+n.ind.cov)) #array cols are Obs ID,Y,Obs species, True species,Obs covariates,Distance,Ind covariates
	G.transect=matrix(0,n.species,n.transects)  #number of groups by transect; each row gives results for separate species
	n.Records=G.transect #number of records by transect (=G.transect*n.Observers)
	N.transect=G.transect #total abundance by transect

	for(isp in 1:n.species){
		for(itrans in 1:n.transects){
			cur.gt0=sum(Dat.num[,"Transect"]==itrans & Dat.num[,"True.sp"]==isp)
			if(cur.gt0>0){
				Cur.dat=Dat.num[which(Dat.num[,"Transect"]==itrans & Dat.num[,"True.sp"]==isp),3:ncol(Dat.num)]
				Data[isp,itrans,1:nrow(Cur.dat),1:(3+n.obs.cov)]=as.matrix(Cur.dat[,1:(3+n.obs.cov)])
				Data[isp,itrans,1:nrow(Cur.dat),5+n.obs.cov]=as.matrix(Cur.dat[,"Distance"])
				#fill distances for unobserved
			  if((M[isp,itrans]-nrow(Cur.dat))<=0)cat(paste("Error: M dimension for species ",isp," transect ",itrans," too small; Increase M \n"))
				if(i.binned==1)Data[isp,itrans,(nrow(Cur.dat)+1):M[isp,itrans],5+n.obs.cov]=rep(sample(c(1:n.bins),size=(M[isp,itrans]-nrow(Cur.dat))/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
				else Data[isp,itrans,(nrow(Cur.dat)+1):M[isp,itrans],5+n.obs.cov]=rep(runif((M[isp,itrans]-nrow(Cur.dat))/n.Observers[itrans]),each=n.Observers[itrans])
				#fill individual covariate values for (potential) animals that weren't observed
				if(n.ind.cov>0){
					Data[isp,itrans,1:nrow(Cur.dat),(6+n.obs.cov):(6+n.obs.cov+n.ind.cov-1)]=as.matrix(Cur.dat[,(6+n.obs.cov):(6+n.obs.cov+n.ind.cov-1)])
					for(icov in 1:n.ind.cov){
						rsamp=switch_sample(n=(M[isp,itrans]-nrow(Cur.dat))/n.Observers[itrans],pdf=Cov.prior.pdf[isp,icov],cur.par=Cov.prior.parms[isp,1:Cov.prior.n[isp,icov],icov],RE=0) 
						Data[isp,itrans,(nrow(Cur.dat)+1):M[isp,itrans],5+n.obs.cov+icov]=rep(rsamp,each=n.Observers[itrans])
					}
				}
				n.Records[isp,itrans]=nrow(Cur.dat)
				G.transect[isp,itrans]=nrow(Cur.dat)/n.Observers[itrans]		#initialize abundance in each transect area to be = to total number of animals observed
				#fill species for latent animals				
				Data[isp,itrans,(n.Records[isp,itrans]+1):max.M,4]=rep(sample(c(1:n.species),(max.M-n.Records[isp,itrans])/n.Observers[itrans],replace=TRUE),each=n.Observers[itrans])
				if(grps==FALSE)N.transect[isp,itrans]=G.transect[isp,itrans]
				else N.transect[isp,itrans]=sum(Cur.dat[,"Group"])/n.Observers[itrans]
			}
			else{
				if(i.binned==1)Data[isp,itrans,1:M[isp,itrans],5+n.obs.cov]=rep(sample(c(1:n.bins),size=M[isp,itrans]/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
				else Data[isp,itrans,1:M[isp,itrans],5+n.obs.cov]=rep(runif(M[isp,itrans]/n.Observers[itrans]),each=n.Observers[itrans])
				if(n.ind.cov>0){
					for(icov in 1:n.ind.cov){
						rsamp=switch_sample(n=M[isp,itrans]/n.Observers[itrans],pdf=Cov.prior.pdf[isp,icov],cur.par=Cov.prior.parms[isp,1:Cov.prior.n[isp,icov],icov],RE=0)
						Data[isp,itrans,1:M[isp,itrans],5+n.obs.cov+icov]=rep(rsamp,each=n.Observers[itrans])
					}
				}
				#fill species
				Data[isp,itrans,1:max.M,4]=rep(sample(c(1:n.species),max.M/n.Observers[itrans],replace=TRUE),each=n.Observers[itrans])
				G.transect[isp,itrans]=0
				n.Records[isp,itrans]=0
				N.transect[isp,itrans]=0
			}
			#fill observer ids
			Data[isp,itrans,(n.Records[isp,itrans]+1):max.M,1]=rep(Observers.num[1:n.Observers[itrans],itrans],(max.M-n.Records[isp,itrans])/n.Observers[itrans])
			#fill observer covariates
			if(n.obs.cov>0){
				for(icov in 1:n.obs.cov){
					Data[isp,itrans,1:max.M,4+icov]=rep(Obs.cov[1:n.Observers[itrans],itrans,icov],max.M/n.Observers[itrans])
				}
			}
		}
	}
	#for debugging, when data aug is set to truth (keep 0's when simulating data), need to put unobserved animals after observed animals
#	if(DEBUG==TRUE){
#		for(isp in 1:n.species){
#			for(itrans in 1:n.transects){
#				if(G.transect[isp,itrans]>0){
#					Cur.dat=Data[isp,itrans,1:n.Records[isp,itrans],]
#					if(n.Observers[itrans]==2){
#						Obs.ind=matrix(Cur.dat[,2],2,G.transect[isp,itrans])
#						Obs.ind[1,]=apply(Obs.ind,2,'max')
#						Obs.ind[2,]=Obs.ind[1,]
#						Obs.ind=as.vector(Obs.ind)
#					}
#					else Obs.ind=Cur.dat[,2]
#					Data[isp,itrans,1:n.Records[isp,itrans],]=rbind(Cur.dat[which(Obs.ind==1),],Cur.dat[which(Obs.ind==0),])			
#				}
#			}
#		}
#	}
#	
	#set starting 'observed' species to zero for nondetections
    for(isp in 1:n.species){
		for(itrans in 1:n.transects){
			Data[isp,itrans,,3]=Data[isp,itrans,,2]*Data[isp,itrans,,3]
		}
	}
	
	#set all 'true species' values = to defining entry in species array
 	for(isp in 1:n.species){
		Data[isp,,,4]=isp
	}

	stacked.names=c(colnames(Dat)[3:4],"Obs.species","Species",colnames(Dat)[6:ncol(Dat)])
	Stacked=stack_data(Data[1,,,],G.transect[1,]*n.Observers,n.transects,stacked.names,factor.ind) #a stacked form of detection data for updating beta parameters
	if(n.species>1)for(isp in 2:n.species)Stacked=rbind(Stacked,stack_data(Data[isp,,,],G.transect[isp,]*n.Observers,n.transects,stacked.names,factor.ind))

	#determine levels for each factor variable to help in assembling compatible DMs for smaller datasets 
	# (stored in list object named 'Levels')
  factor.ind["Species"]=TRUE
	factor.cols=which(factor.ind[stacked.names]==TRUE) 
	if(length(factor.cols)>0 & is.null(Levels)[1]==1){
		Temp=Stacked[,factor.cols]
		Levels=eval(parse(text=paste('list(',colnames(Temp)[1],'=levels(Temp[,1]))',sep='')))
		if(length(factor.cols)>1){
			for(icol in 2:length(factor.cols)){
				eval(parse(text=paste('Levels$',colnames(Temp)[icol],'=levels(Temp[,icol])',sep='')))	
			}
		}		
	}

	N.hab.pois.par=rep(0,n.species)
	DM.hab.pois=vector('list',n.species)
	if(1==1){
		if(is.null(Hab.cov)|Hab.pois.formula[[1]]==~1){
			DM.hab.pois[[1]]=as.matrix(rep(1,S),ncol=1)
			colnames(DM.hab.pois[[1]])="Intercept"
		}
		else DM.hab.pois[[1]]=model.matrix(Hab.pois.formula[[1]],data=Hab.cov)
	}
	N.hab.pois.par[1]=ncol(DM.hab.pois[[1]])
	if(n.species>1){
		for(i in 2:n.species){  #create design matrices for each species. e.g., name for first species will be DM.hab.pois1
			if(is.null(Hab.cov)|Hab.pois.formula[[i]]==~1){
				DM.hab.pois[[i]]=as.matrix(rep(1,S),ncol=1)
				colnames(DM.hab.pois[[i]])="Intercept"
			}
			else DM.hab.pois[[i]]=model.matrix(Hab.pois.formula[[i]],data=Hab.cov)
			N.hab.pois.par[i]=ncol(DM.hab.pois[[i]])
		}
	}
  DM.hab.bern=NULL
  N.hab.bern.par=NULL
	if(ZIP){
	  DM.hab.bern=vector('list',n.species)
	  N.hab.bern.par=rep(0,n.species)
	  if(1==1){
	    if(is.null(Hab.cov)|Hab.bern.formula[[1]]==~1){
	      DM.hab.bern[[1]]=as.matrix(rep(1,S),ncol=1)
	      colnames(DM.hab.bern[[1]])="Intercept"
	    }
	    else DM.hab.bern[[1]]=model.matrix(Hab.bern.formula[[1]],data=Hab.cov)
	  }
	  N.hab.bern.par[1]=ncol(DM.hab.bern[[1]])
	  if(n.species>1){
	    for(i in 2:n.species){  #create design matrices for each species. e.g., name for first species will be DM.hab.bern1
	      if(is.null(Hab.cov)|Hab.bern.formula[[i]]==~1){
	        DM.hab.bern[[i]]=as.matrix(rep(1,S),ncol=1)
	        colnames(DM.hab.bern[[i]])="Intercept"
	      }
	      else DM.hab.bern[[i]]=model.matrix(Hab.bern.formula[[i]],data=Hab.cov)
	      N.hab.bern.par[i]=ncol(DM.hab.bern[[i]])
	    }
	  }
	}	  
  
  DM.det=get_mod_matrix(Cur.dat=Stacked,stacked.names,factor.ind,Det.formula,Levels)
	
	#now, deal with misID parameters
	N.par.misID=NULL
	if(misID==TRUE){
		N.par.misID=rep(0,max(misID.mat))
		for(i in 1:max(misID.mat)){
			N.par.misID[i]=ncol(model.matrix(misID.models[[i]],data=as.data.frame(Stacked)))
		}
	}

	Par=generate_inits_misID(DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,DM.det=DM.det,N.hab.pois.par=N.hab.pois.par,N.hab.bern.par=N.hab.bern.par,G.transect=G.transect,Area.trans=Area.trans,Area.hab=Area.hab,Mapping=Mapping,point.ind=point.ind,spat.ind=spat.ind,grp.mean=Cov.prior.parms[,1,1],misID=misID,misID.mat=misID.mat,N.par.misID=N.par.misID)	
	if(is.null(Inits)==FALSE){  #replace random inits with user provided inits for all parameters specified
		I.init=names(Inits)
		for(ipar in 1:length(I.init)){
			eval(parse(text=paste("Par$",names(Inits)[ipar],"=Inits$",names(Inits[ipar]))))
		}
	}
  if(point.ind==TRUE)Par$cor=0
	#start Nu out at a compatible level
    for(isp in 1:n.species)Par$Nu[isp,]=DM.hab.pois[[isp]]%*%Par$hab.pois[isp,1:N.hab.pois.par[isp]]
  if(length(Par$tau.nu)!=n.species)cat('Error: length of initial value vector for tau.nu should be equal to # of species')
	
	#get initial individual covariate parameter values
	Par$Cov.par=Cov.prior.parms 
	for(i in 1:n.ind.cov){	
		for(j in 1:n.species){
			if(Cov.prior.fixed[j,i]==1)Par$Cov.par[j,,i]=Cov.prior.parms[j,,i]
			else{
				temp=switch_sample_prior(Cov.prior.pdf[j,i],Cov.prior.parms[j,,i])
				Par$Cov.par[j,1:length(temp),i]=temp
			}
		}
	}
	
	dist.pl=5+n.obs.cov
	
	n.hab.cov=ifelse(is.null(Hab.cov)==1 | length(Hab.cov)==1,0,ncol(Hab.cov))
		
	i.Covered=c(1:S)%in%Mapping
	Covered.area=rep(0,S)
	for(i in 1:S){
		if(i.Covered[i]==1){
			Covered.area[i]=sum(Area.trans[which(Mapping==i)])
		}
	}
	
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q=Matrix(Q)	

	Meta=list(n.transects=n.transects,n.species=n.species,S=S,spat.ind=spat.ind,Area.hab=Area.hab,Area.trans=Area.trans,
			Adj=Adj,Mapping=Mapping,Covered.area=Covered.area,n.Observers=n.Observers,M=M,stacked.names=stacked.names,
			factor.ind=factor.ind,Det.formula=Det.formula,detect=detect,Levels=Levels,i.binned=i.binned,dist.pl=dist.pl,
			G.transect=G.transect,N.transect=N.transect,grps=grps,n.bins=n.bins,Bin.length=Bin.length,n.ind.cov=n.ind.cov,
			Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,point.ind=point.ind,last.ind=last.ind,cor.const=cor.const,fix.tau.nu=fix.tau.nu,
			srr=srr,srr.tol=srr.tol,misID=misID,misID.models=misID.models,misID.mat=misID.mat,misID.symm=misID.symm,N.par.misID=N.par.misID,N.hab.pois.par=N.hab.pois.par,N.hab.bern.par=N.hab.bern.par,post.loss=post.loss)
	
	if(adapt==TRUE){
		cat('\n Beginning adapt phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$adapt,adapt=1,Control=Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Out$Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	else{
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	Out	
}
