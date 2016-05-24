#' Function for MCMC analysis 
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		"det": a vector giving the current iteration's linear model parameters for the detection model;
#' 		"hab.pois": a vector giving the current iteration's linear model parameters for Poisson abundance intensity; each row gives parameters for a particular species
#'   	"hab.bern": a vector giving the current iteration's linear model parameters for Bernoulli part of ZIP model for abundance (if Meta$ZIP=TRUE)
#' 		"cor": a correlation parameter for detections that's an increasing function of distance (correlation at the maximum distance);
#' 		"Nu": a vector giving the log of the abundance intensity for each strata;
#'    "Eta.pois": If Meta$spat.ind==FALSE, spatial random effects for Poisson abundance model; one for each cell and for each species
#'    "Eta.bern": If Meta$spat.ind==FALSE & Meta$ZIP=TRUE, spatial random effects for Bernoulli abundance model; one for each cell and for each species
#'	  "tau.eta.pois": If Meta$spat.ind==FALSE, precision for spatial ICAR model(s) for the Poisson component
#'    "tau.eta.bern": If Meta$spat.ind==FALSE & Meta$ZIP=TRUE, precision for spatial ICAR model(s) for the Bernoulli component
#'	  "tau.nu": Precision for Nu (overdispersion relative to the Poisson distribution)
#' 		"G": a vector giving the number of groups of animals in each strata; 
#' 		"N": a vector giving the number of animals in each strata
#' 		"MisID": a list, each entry i of which is a vector holding the parameters for the ith misID model
#' 		"Cov.par": an (n.species X n X n.ind.cov)  array holding parameters of individual covariate distributions.
#' @param Data   A four dimensional array; the first dimension gives species, the second gives the transect, the third indexes a (possible) observation, 
#' 			and the fourth dimension gives observations and covariates associated with a given animal.
#' 			These final columns are: Observer ID,Y(observation=0/1),Observed species,Obs covariates,Distance,Ind covariates
#' @param cur.iter   Number of iterations to run
#' @param adapt	If adapt==TRUE, run MCMC in adapt mode, optimizing MCMC proposal distributions prior to primary MCMC
#' @param Control	A list object including the following objects:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal 
#' 				ranges prior to final MCMC run; 
#'	"MH.cor": Metropolis-hastings tuning parameter for updating the correlation parameter (if Meta$point.ind==TRUE);
#'	"MH.nu": MH tuning parameter for Nu parameters (Langevin-Hastings multivariate update);
#'	"RJ.N": A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect
#'  "iter.fix.N"  Number of iterations to skip RJMCMC step 
#' @param DM.hab.pois	A design matrix for the Poisson model for abundance intensity (log scale)
#' @param DM.hab.bern If Meta$ZIP=TRUE, a design matrix for the Bernoulli zero model (probit scale)
#' @param DM.det	A design matrix for the probit of detection probability
#' @param Q			An inverse precision matrix for the spatial ICAR process
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following objects
#'	"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'  "b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'	"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
#'	"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
#'	"beta.tau": Prior precision for regression coefficients 
#'  "misID.mu": a list vector, each entry gives normal prior means for misID regression coefficients for the corresponding model in Meta$misID.mat (can be set to null if no misID)
#'  "misID.sd": a list vector, each entry gives normal prior sd for misID regression coefficients for the corresponding model in Meta$misID.mat (can be set to null if no misID)
#' @param Meta	A list object giving a number of other features of the dataset, including:
#' 	"n.transects"	Number of transects
#'  "n.species"     Number of species
#' 	"S"				Number of strata cells
#'  "spat.ind"		Indicator for spatial dependence
#'  "Area.hab"		Vector giving relative area covered by each strata
#'  "Area.trans"	Vector giving fraction of area of relevant strata covered by each transect
#' 	"Adj"			Adjacency matrix giving connectivity of spatial grid cells
#'  "Mapping" 		Vector mapping each transect into a parent strata
#'  "Covered.area"	Vector giving the fraction of each strata covered by transects
#' 	"n.Observers"	Vector giving the number of observers that operated on each transect
#'  "M"   Matrix with species-specific rows giving maximum possible value for number of groups present in each transect (in practice just set high enough that values at M and above are never sampled during MCMC) and can be fine tuned as needed#'  "stacked.names" Character vector giving column names for the dataset
#'  "factor.ind"	Indicator vector specifying whether data columns are factors (1) or continuous (0)
#'  "detect"		If TRUE, detection parameters are estimated; if FALSE assumes a census
#'  "Det.formula"	a formula object specifying the model for the detection process
#'  "Levels"		a list object, whose elements are comprised of detection model names; each element gives total # of levels in the combined dataset
#'  "i.binned"		indicator for whether distances are recorded in bins (1) or are continuous (0)
#'  "dist.pl"		gives the column in Data where distances are located	
#'  "G.transect"	vector holding current number of groups of animals present in area covered by each transect		
#'  "N.transect"    vector holding current number of animals present in covered area by each transect
#'  "grps"			indicator for whether observations are for groups rather than individuals
#'  "n.bins"		number of distance bins (provided i.binned=1)
#'  "Bin.length"	vector giving relative size of distance bins 
#'  "n.ind.cov" 	Number of individual covariates (distance is not included in this total, but group size is)
#'  "Cov.prior.pdf" character vector giving the probability density function associated with each individual covariate (type ? hierarchical_DS for more info)
#'  "Cov.prior.parms"	An (n.species X n X n.ind.cov) array providing "pseudo-prior" parameters for individual covarate distributions (only the first row used if a signle parameter distribution)
#'  "Cov.prior.fixed" indicator vector for whether parameters of each covariate distribution should be fixed within estimation routine
#'  "Cov.prior.n" 	(#species X #covariates) Matrix giving number of parameters in each covariate pdf 
#'  "ZIP"  If TRUE, fit a ZIP model to abundance
#'  "point.ind"		Indicator for whether point independence assumed (if no, then no correlation modeled b/w multiple observers as function of distance)
#'  "last.ind" If TRUE (and point.ind=TRUE), point independence operates by assuming 0 dependence at the farthest bin
#'  "cor.const" If TRUE, forces estimates of correlation associated with point independence to be positive if last.ind==FALSE or negative if last.ind==TRUE (default is FALSE)
#'  "fix.tau.nu"	Indicator for whether tau.nu should be fixed (1) or estimated(0)
#'  "srr"			Indicator for whether a spatially restricted regression model should be employed (1) or not (0)
#'  "srr.tol"		Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation
#'  "misID"			If TRUE, misidentification of species is modeled
#'  "misID.mat"     With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#'  "misID.models"  A formula vector providing linar model-type formulas for each positive value of misID.mat. 
#'  "misID.symm"    If TRUE, classification probabilities assumed to be symmetric (e.g. pi^{2|1}=pi^{1|2})
#'  "N.par.misID"   A vector specifying the number of parameters needed for each misID model
#'  "N.hab.pois.par"	    A vector specifying the number of parameters needed for each species' Poisson abundance model
#'  "N.hab.bern.par"    If fitting a ZIP model, this vector specifying the number of parameters needed for each species' Bernoulli zero model
#'  "post.loss"  If TRUE, observed and predicted detections are compiled for posterior predictive loss 
#' @return returns a list with the following objects: 
#' 	"MCMC": An 'mcmc' object (see 'coda' R package) containing posterior samples;
#'  "Accept": A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  "Control": A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used) 
#'  "Obs.N":  Records latent abundance in each transect; dimension is (n.species X # samples X # transects)
#'  "Pred.N": Posterior predictive distribution for abundance in each transect; obtained by sampling a Poisson distribution given current parameter values (with possible zero inflation)
#'  "Post": Holds posterior samples for strata specific group sizes ("Post$G") and abundance ("Post$N")
#'  "Obs.det":  if Meta$post.loss=TRUE, an array holding observed detection types for posterior predictive loss calculations dim = c(n.transects,n.obs.types,n.obs.types) 
#'  "Pred.det": if Meta$post.loss=TRUE, an array holding predicted detection types for posterior predictive loss calculations dim = c(n.mcmc.iter,n.transects,n.obs.types,n.obs.types)
#' @export
#' @import Matrix
#' @keywords areal, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn

mcmc_ds<-function(Par,Data,cur.iter,adapt,Control,DM.hab.pois,DM.hab.bern=NULL,DM.det,Q,Prior.pars,Meta){	
	#require(mvtnorm)
	#require(Matrix)
	#require(truncnorm)
  SMALL=10^{-20}
	Lam.index=c(1:Meta$S)
	if(Meta$i.binned==0)dist.mult=1
	if(Meta$i.binned==1)dist.mult=1/(Meta$n.bins-1)
	n.beta.det=ncol(DM.det)
	n.Records=t(t(Meta$G.transect)*Meta$n.Observers)
	grp.pl=NULL
	if(Meta$grps==TRUE)grp.pl=which(Meta$stacked.names=="Group")
	
	#initialize G.obs (number of groups observed per transect)
	G.obs=Meta$G.transect  #number of groups observed by at least one observer
	for(isp in 1:Meta$n.species){
		for(itrans in 1:Meta$n.transects){
			Tmp=matrix(Data[isp,itrans,1:n.Records[isp,itrans],2],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
			G.obs[isp,itrans]=sum(apply(Tmp,1,'sum')>0)
		}
	}
	g.tot.obs=colSums(G.obs)%*%Meta$n.Observers  #total number of observations of animals seen at least once
	n.obs.cov=Meta$dist.pl-4
	
	n.samp.misID=max(1,round(0.05*sum(G.obs)))  #currently only updating species for 1/10 of population at each iteration
	
	if(Meta$detect==FALSE){
		Meta$Det.formula=~1
		Par$det=10
	}
	
	#initialize Y.tilde (temporarily pretending correlation between observations is zero)
	if(Meta$detect){
		Y.tilde=array(0,dim=c(Meta$n.species,max(Meta$M),Meta$n.transects))
		for(isp in 1:Meta$n.species){
			for(itrans in 1:Meta$n.transects){
        #cat(paste("isp ",isp," itrans ",itrans))
				X=get_mod_matrix(Cur.dat=Data[isp,itrans,,],Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
				ExpY=X%*%Par$det
				Y.tilde[isp,,itrans]=rtruncnorm(max(Meta$M), a=ifelse(Data[isp,itrans,,2]==0,-Inf,0), b=ifelse(Data[isp,itrans,,2]==0,0,Inf), ExpY, 1)		
			}
		}
	}
 
	Z=matrix(1,Meta$n.species,Meta$S)  #need to define for non-ZIP models
  #in case of ZIP model, initialize Z, Z.tilde
  if(Meta$ZIP){
    #start all zeros as arising from Bernoulli component
    Z[Par$G==0]=0    
    Z.tilde=matrix(0,Meta$n.species,Meta$S)
    G.gt0=(Par$G>0)
    G.eq0=(Par$G==0)
    for(isp in 1:Meta$n.species){
      n.gt0=sum(G.gt0[isp,])
      n.eq0=sum(G.eq0[isp,])
      ExpZ=DM.hab.bern[[isp]]%*%Par$hab.bern[isp,]
      if(n.gt0>0)Z.tilde[isp,which(G.gt0[isp,]==1)]=rtruncnorm(n.gt0,a=0,b=Inf,ExpZ[G.gt0[isp,]],1)
      if(n.eq0>0)Z.tilde[isp,which(G.eq0[isp,]==1)]=rtruncnorm(n.eq0,a=-Inf,b=0,ExpZ[G.eq0[isp,]],1)
    }
  }

	#initialize lambda
	Lambda=matrix(0,Meta$n.species,Meta$S)
	Lambda.trans=matrix(0,Meta$n.species,Meta$n.transects)
	for(isp in 1:Meta$n.species){
		Lambda[isp,]=exp(Par$Nu[isp,])*Meta$Area.hab
		Lambda.trans[isp,]=Lambda[isp,Meta$Mapping]*Meta$Area.trans
	}
	grp.lam=rep(0,Meta$n.species)
	
	
	#initialize statistics/matrices needed for MCMC updates
	XpXinv.pois=vector('list',Meta$n.species)
	XpXinvXp.pois=XpXinv.pois
	for(isp in 1:Meta$n.species){
		XpXinv.pois[[isp]]=solve(crossprod(DM.hab.pois[[isp]]))
		XpXinvXp.pois[[isp]]=XpXinv.pois[[isp]]%*%t(DM.hab.pois[[isp]])
	}
  
  if(Meta$ZIP){
    XpXinv.bern=vector('list',Meta$n.species)
    XpXinvXp.bern=XpXinv.bern
    for(isp in 1:Meta$n.species){
      XpXinv.bern[[isp]]=solve(crossprod(DM.hab.bern[[isp]]))
      XpXinvXp.bern[[isp]]=XpXinv.bern[[isp]]%*%t(DM.hab.bern[[isp]])
    }
  }

	if(Meta$srr){
		L.t.pois=XpXinv.pois
		L.pois=L.t.pois
		Qt.pois=L.t.pois
		cross.L.pois=L.t.pois
		Theta.pois=L.t.pois
		N.theta.pois=rep(0,Meta$n.species)
		for(isp in 1:Meta$n.species){
			P.c=diag(Meta$S)-DM.hab.pois[[isp]]%*%solve(crossprod(DM.hab.pois[[isp]]),t(DM.hab.pois[[isp]]))
			Omega=(P.c%*%Meta$Adj%*%P.c)*(Meta$S/sum(Meta$Adj))
			Eigen=eigen(Omega)
			if(max(Eigen$values)<Meta$srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
			Ind=which(Eigen$values>Meta$srr.tol)
			L.t.pois[[isp]]=Eigen$vectors[,Ind]
			cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
			L.pois[[isp]]=t(L.t.pois[[isp]])
			Qt.pois[[isp]]=L.pois[[isp]]%*%Q%*%L.t.pois[[isp]]
			cross.L.pois[[isp]]=L.pois[[isp]]%*%L.t.pois[[isp]]
			N.theta.pois[isp]=nrow(Qt.pois[[isp]])
			Theta.pois[[isp]]=rnorm(N.theta.pois[isp],0,sqrt(1/Par$tau.eta.pois[isp]))
		}
    if(Meta$ZIP){
      L.t.bern=XpXinv.bern
      L.bern=L.t.bern
      Qt.bern=L.t.bern
      cross.L.bern=L.t.bern
      Theta.bern=L.t.bern
      N.theta.bern=rep(0,Meta$n.species)
      for(isp in 1:Meta$n.species){
        P.c=diag(Meta$S)-DM.hab.bern[[isp]]%*%solve(crossprod(DM.hab.bern[[isp]]),t(DM.hab.bern[[isp]]))
        Omega=(P.c%*%Meta$Adj%*%P.c)*(Meta$S/sum(Meta$Adj))
        Eigen=eigen(Omega)
        if(max(Eigen$values)<Meta$srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
        Ind=which(Eigen$values>Meta$srr.tol)
        L.t.bern[[isp]]=Eigen$vectors[,Ind]
        cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
        L.bern[[isp]]=t(L.t.bern[[isp]])
        Qt.bern[[isp]]=L.bern[[isp]]%*%Q%*%L.t.bern[[isp]]
        cross.L.bern[[isp]]=L.bern[[isp]]%*%L.t.bern[[isp]]
        N.theta.bern[isp]=nrow(Qt.bern[[isp]])
        Theta.bern[[isp]]=rnorm(N.theta.bern[isp],0,sqrt(1/Par$tau.eta.bern[isp]))
      }
    }
	}
  
  
	Sampled=unique(Meta$Mapping)
	n.unique=length(Sampled)
	Sampled.area.by.strata=rep(0,n.unique)
	for(i in 1:Meta$n.transects)Sampled.area.by.strata[Sampled==Meta$Mapping[i]]=Sampled.area.by.strata[which(Sampled==Meta$Mapping[i])]+Meta$Area.trans[i]
	
	#initialize MCMC, Acceptance rate matrices
	mcmc.length=(Control$iter-Control$burnin)/Control$thin
  MCMC=list(N.tot=matrix(0,Meta$n.species,mcmc.length),N=array(0,dim=c(Meta$n.species,mcmc.length,Meta$S)),G=array(0,dim=c(Meta$n.species,mcmc.length,Meta$S)),Hab.pois=array(0,dim=c(Meta$n.species,mcmc.length,ncol(Par$hab.pois))),Det=data.frame(matrix(0,mcmc.length,length(Par$det))),cor=rep(0,mcmc.length),tau.eta.pois=matrix(0,Meta$n.species,mcmc.length),tau.nu=matrix(0,Meta$n.species,mcmc.length),Cov.par=array(0,dim=c(Meta$n.species,mcmc.length,length(Par$Cov.par[1,,]))))
  if(Meta$ZIP){
    MCMC$Hab.bern=array(0,dim=c(Meta$n.species,mcmc.length,ncol(Par$hab.bern)))
    MCMC$tau.eta.bern=matrix(0,Meta$n.species,mcmc.length)
  }
  if(Meta$misID){
    MCMC$MisID=vector("list",length(Meta$N.par.misID))
    for(ipar in 1:length(Meta$N.par.misID))MCMC$MisID[[ipar]]=matrix(0,Meta$N.par.misID[ipar],mcmc.length)
	}
  #colnames(MCMC$Hab)=colnames(DM.hab)
	if(Meta$detect)colnames(MCMC$Det)=colnames(DM.det)
	if(Meta$misID){
		Accept=list(cor=0,N=matrix(0,Meta$n.species,Meta$n.transects),Nu=matrix(0,Meta$n.species,n.unique),MisID=vector("list",length(Meta$N.par.misID)))
		for(ipar in 1:length(Meta$N.par.misID))Accept$MisID[[ipar]]=rep(0,Meta$N.par.misID[ipar])
	}
	if(Meta$misID==FALSE)Accept=list(cor=0,N=matrix(0,Meta$n.species,Meta$n.transects),Nu=matrix(0,Meta$n.species,n.unique))
	Pred.N=array(0,dim=c(Meta$n.species,mcmc.length,Meta$n.transects))
	Obs.N=Pred.N

	#check that intial misID parameters are plausible, set multinomial arrays (needed for MH updating)
	n.obs.types=ncol(Meta$misID.mat)
	if(Meta$misID){
		for(isp in 1:Meta$n.species){
			Cur.dat=matrix(0,sum(G.obs[isp,]*Meta$n.Observers),dim(Data)[4])
			ipl=1
			for(itrans in 1:Meta$n.transects){
				if(G.obs[isp,itrans]>0)Cur.dat[ipl:(ipl+G.obs[isp,itrans]*Meta$n.Observers[itrans]-1),]=Data[isp,itrans,1:(G.obs[isp,itrans]*Meta$n.Observers[itrans]),]
				ipl=ipl+G.obs[isp,itrans]*Meta$n.Observers[itrans]
			}
			if(Meta$detect==1)Cur.dat=Cur.dat[-which(Cur.dat[,3]==0),]
			DM=vector('list',n.obs.types)
			XBeta=matrix(0,nrow(Cur.dat),n.obs.types)	
			for(ipl in 1:n.obs.types){  #set up initial parameter values, cell probabilities
				if(Meta$misID.mat[isp,ipl]>0){
					DM[[ipl]]=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[Meta$misID.mat[isp,ipl]]],Levels=Meta$Levels)
					XBeta[,ipl]=exp(DM[[ipl]]%*%Par$MisID[[Meta$misID.mat[isp,ipl]]])
				}
			}
			if(sum(XBeta[,isp]<0 | XBeta[,isp]<apply(XBeta[,-isp],1,'max'))!=0){ 
				cat('\n Error: Initial misidentifiation parameters are not consistent with assumption that correct ID probability is > misID probability for all animals \n')#implies XB bigger than other ID possibilities for all animals
			}
		}
		Multinom.resp=array(0,dim=c(Meta$n.transects,Meta$n.species,n.obs.types))
		Multinom.index=matrix(0,Meta$n.transects,Meta$n.species) 
		for(itrans in 1:Meta$n.transects){
			for(isp in 1:Meta$n.species){
			  if(Meta$G.transect[isp,itrans]>0){
				  for(iobs in 1:(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans])){
					  if(Data[isp,itrans,iobs,3]>0)Multinom.resp[itrans,isp,Data[isp,itrans,iobs,3]]=Multinom.resp[itrans,isp,Data[isp,itrans,iobs,3]]+1	
				  }
			  }
			Multinom.index[itrans,]=rowSums(Multinom.resp[itrans,,])
			}
		}
	}
	Obs.det=NA
  Pred.det=NA
	if(Meta$post.loss){ #calculate observed counts of different detection types, initialize prediction arrays
    Obs.det=array(0,dim=c(Meta$n.transects,Meta$n.species+2,Meta$n.species+2)) #row/col=1 is 'undetected'
    Pred.det=array(0,dim=c(mcmc.length,Meta$n.transects,Meta$n.species+2,Meta$n.species+2))
    for(itrans in 1:Meta$n.transects){
      if(Meta$n.Observers[itrans]==1){ #in this case, just fill first column
        for(isp in 1:Meta$n.species){
          if(G.obs[isp,itrans]>0){
            for(iind in 1:G.obs[isp,itrans]){
              Obs.det[itrans,Data[isp,itrans,iind,3]+1,1]=Obs.det[itrans,Data[isp,itrans,iind,3]+1,1]+1
            }
          }
        }
      }
      else{  
        for(isp in 1:Meta$n.species){
          if(G.obs[isp,itrans]>0){
            for(iind in 1:G.obs[isp,itrans]){
              Obs.det[itrans,Data[isp,itrans,iind*Meta$n.Observers[itrans]-1,3]+1,Data[isp,itrans,iind*Meta$n.Observers[itrans],3]+1]=Obs.det[itrans,Data[isp,itrans,iind*Meta$n.Observers[itrans]-1,3]+1,Data[isp,itrans,iind*Meta$n.Observers[itrans],3]+1]+1
            }
          }
        }
      }      
    }
  }
	
	#initialize random effect matrices for individual covariates if required
	if(sum(1-Meta$Cov.prior.fixed)>0)RE.cov=array(0,dim=c(Meta$n.species,Meta$n.transects,max(Meta$M),Meta$n.ind.cov))
	
	PROFILE=FALSE  #outputs time associated with updating different groups of parameters
	DEBUG=FALSE
	if(DEBUG){
		Par$misID[[1]]=2
		Par$misID[[2]]=1
		Par$misID[[3]]=3
	}
	st <- Sys.time()
	##################################################
	############   Begin MCMC Algorithm ##############
	##################################################
	for(iiter in 1:cur.iter){
		#cat(paste('\n ', iiter))
		for(isp in 1:Meta$n.species){		

		  ########## update abundance parameters at the strata scale   ################
		  Hab.pois=Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]
		  Eta.pois=Par$Eta.pois[isp,]
		  
      if(Meta$ZIP){
		    Hab.bern=Par$hab.bern[isp,1:Meta$N.hab.bern.par[isp]]
		    Eta.bern=Par$Eta.bern[isp,]
      }
        
		  #update Z,Z.tilde if ZIP model specified
		  if(Meta$ZIP){
		    Z[isp,which(Par$G[isp,]>0)]=1
		    Which.0=which(Par$G[isp,]==0)
		    Cur.p=pnorm(DM.hab.bern[[isp]]%*%Hab.bern+Eta.bern) #bernoulli success prob
		    Cur.pois0=Cur.p*exp(-Meta$Area.hab*exp(Par$Nu[isp,]))
		    Cur.p=Cur.pois0/((1-Cur.p)+Cur.pois0)
		    if(length(Which.0)>0)Z[isp,which(Par$G[isp,]==0)]=rbern(length(Which.0),Cur.p[Which.0])
		    Z.tilde[isp,]=rtruncnorm(Meta$S,a=ifelse(Z[isp,]==1,0,-Inf),b=ifelse(Z[isp,]==1,Inf,0),DM.hab.bern[[isp]]%*%Hab.bern+Eta.bern)       
		  }
		  
		  #update nu parameters (log lambda)
	 	  #1) for sampled cells for which z-tilde>0 (if ZIP)
		  Mu=DM.hab.pois[[isp]]%*%Hab.pois+Eta.pois
		  G.sampled=rep(0,n.unique) #total number of groups currently in each sampled strata
		  for(i in 1:Meta$n.transects)G.sampled[Sampled==Meta$Mapping[i]]=G.sampled[Sampled==Meta$Mapping[i]]+Meta$G.transect[isp,i]
      for(i in 1:n.unique){
        if(!Meta$ZIP|Z[isp,Sampled[i]]==1){
		      prop=Par$Nu[isp,Sampled[i]]+runif(1,-Control$MH.nu[isp,i],Control$MH.nu[isp,i])				
		      old.post=dnorm(Par$Nu[isp,Sampled[i]],Mu[Sampled[i]],1/sqrt(Par$tau.nu[isp]),log=TRUE)+dpois(G.sampled[i],Sampled.area.by.strata[i]*Meta$Area.hab[Sampled[i]]*exp(Par$Nu[isp,Sampled[i]]),log=TRUE)
		      new.post=dnorm(prop,Mu[Sampled[i]],1/sqrt(Par$tau.nu[isp]),log=TRUE)+dpois(G.sampled[i],Sampled.area.by.strata[i]*Meta$Area.hab[Sampled[i]]*exp(prop),log=TRUE)
		      if(runif(1)<exp(new.post-old.post)){
		        Par$Nu[isp,Sampled[i]]=prop
		        Accept$Nu[isp,i]=Accept$Nu[isp,i]+1
		      }
        }
		  }
		
			#2) simulate nu for areas not sampled
			Par$Nu[isp,-Sampled]=rnorm(Meta$S-n.unique,Mu[-Sampled],1/sqrt(Par$tau.nu[isp]))
		  #3) if ZIP, sample nu for areas where z.tilde<0
      if(Meta$ZIP){
        which.Z.eq0=which(Z[isp,]==0)   
        sampled.Z0=which.Z.eq0[which.Z.eq0%in%Sampled]
        if(length(sampled.Z0)>0)Par$Nu[isp,sampled.Z0]=rnorm(length(sampled.Z0),Mu[sampled.Z0],1/sqrt(Par$tau.nu[isp]))
      }
		  if(PROFILE==TRUE){
		    cat(paste("Nu: ", (Sys.time()-st),'\n'))
		    st=Sys.time()
		  } 
			
      #update spatial random effects
			if(Meta$spat.ind==FALSE){
				if(Meta$srr==FALSE){
					#update eta parameters (spatial random effects)
					V.eta.inv <- Par$tau.nu[isp]*diag(Meta$S) + Par$tau.eta.pois[isp]*Q
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*(Par$Nu[isp,]-DM.hab.pois[[isp]]%*%Hab.pois))		
					Par$Eta.pois[isp,]<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
					Par$Eta.pois[isp,]=Par$Eta.pois[isp,]-mean(Par$Eta.pois[isp,])  #centering
					#update tau_eta  (precision of spatial process)
					Par$tau.eta.pois[isp] <- rgamma(1, (Meta$S-1)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta.pois[isp,], Q %*% Par$Eta.pois[isp,])*0.5) + Prior.pars$b.eta)
          if(Meta$ZIP){
            Hab.bern=Par$hab.bern[isp,]
            V.eta.inv <- diag(Meta$S) + Par$tau.eta.bern[isp]*Q
            M.eta <- solve(V.eta.inv, Z.tilde[isp,]-DM.hab.bern[[isp]]%*%Hab.bern)		
            Par$Eta.bern[isp,]<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
            Par$Eta.bern[isp,]=Par$Eta.bern[isp,]-mean(Par$Eta.bern[isp,])  #centering
            #update tau_eta  (precision of spatial process)
            Par$tau.eta.bern[isp] <- rgamma(1, (Meta$S-1)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta.bern[isp,], Q %*% Par$Eta.bern[isp,])*0.5) + Prior.pars$b.eta)            
					}
				}
				else{
					#Update Theta
					Dat.minus.Exp=Par$Nu[isp,]-DM.hab.pois[[isp]]%*%Hab.pois
					V.eta.inv <- cross.L.pois[[isp]]*Par$tau.nu[isp] + Par$tau.eta.pois[isp]*Qt.pois[[isp]]
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*L.pois[[isp]]%*%Dat.minus.Exp)
					Theta.pois <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(N.theta.pois[isp],0,1))
					Par$Eta.pois[isp,]=as.numeric(L.t.pois[[isp]]%*%Theta.pois)					
					#update tau.eta
					Par$tau.eta.pois[isp] <- rgamma(1, N.theta.pois[isp]*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta.pois, Qt.pois[[isp]] %*% Theta.pois)*0.5) + Prior.pars$b.eta)
          if(Meta$ZIP){
            Hab.bern=Par$hab.bern[isp,]
            Dat.minus.Exp=Z.tilde[isp,]-DM.hab.bern[[isp]]%*%Hab.bern
            V.eta.inv <- cross.L.bern[[isp]] + Par$tau.eta.bern[isp]*Qt.bern[[isp]]
            M.eta <- solve(V.eta.inv, L.bern[[isp]]%*%Dat.minus.Exp)
            Theta.bern <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(N.theta.bern[isp],0,1))
            Par$Eta.bern[isp,]=as.numeric(L.t.bern[[isp]]%*%Theta.bern)	
            #update tau.eta
            Par$tau.eta.bern[isp] <- rgamma(1, N.theta.bern[isp]*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta.bern, Qt.bern[[isp]] %*% Theta.bern)*0.5) + Prior.pars$b.eta)
          }
				}
			}
			#update tau_nu	 (precision for Poisson overdispersion)
			Mu=DM.hab.pois[[isp]]%*%Hab.pois+Par$Eta.pois[isp,]
			if(Meta$fix.tau.nu==FALSE){
        Cur.ind=c(Sampled,which(Z[isp,]==1))
        Cur.ind=Cur.ind[duplicated(Cur.ind)]
				Diff=Par$Nu[isp,Cur.ind]-Mu[Cur.ind]
				Par$tau.nu[isp] <- rgamma(1,length(Cur.ind)/2 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.nu)
			}
			if(PROFILE==TRUE){
				cat(paste("Tau nu: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}

			#translate to lambda scale 
      Lambda[isp,]=exp(Par$Nu[isp,])*Meta$Area.hab  
		  if(Meta$ZIP)Lambda[isp,]=Lambda[isp,]*Z[isp,] 
			Lambda.trans[isp,]=Lambda[isp,Meta$Mapping]*Meta$Area.trans

			if(DEBUG==FALSE){
				#update Betas for habitat relationships
			  Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]=rmvnorm(1,XpXinvXp.pois[[isp]]%*%(Par$Nu[isp,]-Par$Eta.pois[isp,]),XpXinv.pois[[isp]]/(Par$tau.nu[isp]+Prior.pars$beta.tau))
        if(Meta$ZIP)Par$hab.bern[isp,1:Meta$N.hab.bern.par[isp]]=rmvnorm(1,XpXinvXp.bern[[isp]]%*%(Z.tilde[isp,]-Par$Eta.bern[isp,]),XpXinv.bern[[isp]]/(1+Prior.pars$beta.tau))
			}
		  
			########## update group abundance at strata level
			Par$G[isp,]=rpois(Meta$S,Lambda[isp,]*(1-Meta$Covered.area))
			grp.lam[isp]=ifelse(Meta$Cov.prior.pdf[isp,1] %in% c("pois1_ln","poisson_ln"),exp(Par$Cov.par[isp,1,1]+(Par$Cov.par[isp,2,1])^2/2),Par$Cov.par[isp,1,1])
			Par$N[isp,]=Par$G[isp,]+rpois(Meta$S,grp.lam[isp]*Par$G[isp,]) #add the Par$G since number in group > 0 
	    for(ipl in 1:length(Meta$Mapping)){
        Par$G[isp,Meta$Mapping[ipl]]=Par$G[isp,Meta$Mapping[ipl]]+Meta$G.transect[isp,ipl]
			  Par$N[isp,Meta$Mapping[ipl]]=Par$N[isp,Meta$Mapping[ipl]]+Meta$N.transect[isp,ipl]
	    }
		
			if(PROFILE==TRUE){
				cat(paste("Habitat vars, etc.: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			########## update abundance, distances, ind. covariates for observed transects using RJMCMC  #############

			if(Meta$detect & iiter>Control$iter.fix.N){
				for(itrans in 1:Meta$n.transects){
					Sample=c(-Control$RJ.N[isp,itrans]:Control$RJ.N[isp,itrans])
					Sample=Sample[-(Control$RJ.N[isp,itrans]+1)] #don't make proposals where pop size stays the same
					a=sample(Sample,1)
					Sigma=matrix(Par$cor,Meta$n.Observers[itrans],Meta$n.Observers[itrans])
					diag(Sigma)=1
					offdiag=which(Sigma!=1)
					
					if(a>0){ # proposal an addition
						if(((Meta$G.transect[isp,itrans]+a)*Meta$n.Observers[itrans])>=Meta$M[isp,itrans])cat(paste('\n Warning: proposed abundance for transect ',itrans,' species ',isp, '> M; consider increasing M value! \n'))
              #{
							#if(adapt==FALSE)cat(paste('\n Warning: proposed abundance for transect ',itrans,' species ',isp, '> M; consider increasing M value! \n'))
							#else{
						#		temp=floor(Meta$M[isp,itrans]*1.25)
						#		if(temp%%2==1)temp=temp+1
						#		Meta$M[isp,itrans]=min(temp,max(Meta$M[isp,]))
						#	}
						#}
						else{
							Cur.dat=Data[isp,itrans,(n.Records[isp,itrans]+1):(n.Records[isp,itrans]+a*Meta$n.Observers[itrans]),]
							if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
							X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$Det.formula,Levels=Meta$Levels)
							ExpY=X%*%Par$det
							P=c(1:a)
							for(i in 1:a){
								if(Meta$last.ind)Sigma[offdiag]=Par$cor*(Meta$i.binned*(Meta$n.bins-Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])*dist.mult+(1-Meta$i.binned)*(1-Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]))
                else Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
								P[i]=max(pmvnorm(upper=rep(0,Meta$n.Observers[itrans]),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma),SMALL)
							}
              tmp.sum=0
							for(i in 1:a)tmp.sum=tmp.sum-log(Meta$G.transect[isp,itrans]-G.obs[isp,itrans]+i)+log(P[i])
							MH.prob=exp(a*log(Lambda.trans[isp,itrans])+tmp.sum)
							if(runif(1)<MH.prob){
								Meta$G.transect[isp,itrans]=Meta$G.transect[isp,itrans]+a
								n.Records[isp,itrans]=Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]
                if(Meta$grps==FALSE)Meta$N.transect[isp,itrans]=Meta$G.transect[isp,itrans]
								else Meta$N.transect[isp,itrans]=sum(Data[isp,itrans,1:n.Records[isp,itrans],Meta$stacked.names=="Group"])/Meta$n.Observers[itrans]
                Accept$N[isp,itrans]=Accept$N[isp,itrans]+1
								#generate Y-tilde values
								Tmp=matrix(ExpY,a,Meta$n.Observers[itrans],byrow=TRUE)
								Y.tilde.temp=rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,1],sd=1)
								if(Meta$n.Observers[itrans]>1){
									Dist=matrix(Cur.dat[,Meta$dist.pl],a,Meta$n.Observers[itrans],byrow=TRUE)
									Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
									Y.tilde.temp=cbind(Y.tilde.temp,rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,2]+Cor*(Y.tilde.temp-Tmp[,1]),sd=sqrt(1-Cor^2)))
								}
								Y.tilde[isp,(n.Records[isp,itrans]-a*Meta$n.Observers[itrans]+1):n.Records[isp,itrans],itrans]=as.vector(t(Y.tilde.temp))
							}	
						}
					}
					else{ #proposal a deletion
						if((Meta$G.transect[isp,itrans]+a)>=G.obs[isp,itrans]){ #can only delete records where an animal isn't observed
							Cur.dat=Data[isp,itrans,n.Records[isp,itrans]:(n.Records[isp,itrans]+a*Meta$n.Observers[itrans]+1),]
							if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
							X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
							ExpY=X%*%Par$det
							P=c(1:-a)
							for(i in 1:-a){
							  if(Meta$last.ind)Sigma[offdiag]=Par$cor*(Meta$i.binned*(Meta$n.bins-Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])*dist.mult+(1-Meta$i.binned)*(1-Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]))
							  else Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
							  P[i]=pmvnorm(upper=rep(0,Meta$n.Observers[itrans]),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
							}
							tmp.sum=0
							for(i in 1:-a){
								tmp.sum=tmp.sum+log(Meta$G.transect[isp,itrans]-G.obs[isp,itrans]-i+1)-log(P[i])
							}
							MH.prob=exp(a*log(Lambda.trans[isp,itrans])+tmp.sum)
							if(runif(1)<MH.prob){
								Meta$G.transect[isp,itrans]=Meta$G.transect[isp,itrans]+a
								n.Records[isp,itrans]=Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]
								Accept$N[isp,itrans]=Accept$N[isp,itrans]+1
								if(Meta$grps==FALSE)Meta$N.transect[isp,itrans]=Meta$G.transect[isp,itrans]
								else{
                  Meta$N.transect[isp,itrans]=0
                  if(Meta$G.transect[isp,itrans]>0)Meta$N.transect[isp,itrans]=sum(Data[isp,itrans,1:n.Records[isp,itrans],Meta$stacked.names=="Group"])/Meta$n.Observers[itrans]
								}
							}												
						}
					}
					
					#update distances, individual covariates for latent animals not currently in the population
					#fill distances for unobserved
					if(Meta$i.binned==1)Data[isp,itrans,(n.Records[isp,itrans]+1):Meta$M[isp,itrans],Meta$dist.pl]=rep(sample(c(1:Meta$n.bins),size=(Meta$M[isp,itrans]-n.Records[isp,itrans])/Meta$n.Observers[itrans],replace=TRUE,prob=Meta$Bin.length),each=Meta$n.Observers[itrans])
					else Data[isp,itrans,(n.Records[isp,itrans]+1):Meta$M[isp,itrans],Meta$dist.pl]=rep(runif((Meta$M[isp,itrans]-n.Records[isp,itrans])/Meta$n.Observers[itrans]),each=Meta$n.Observers[itrans])
					#fill individual covariate values for (potential) animals that weren't observed
					if(Meta$n.ind.cov>0){
						for(icov in 1:Meta$n.ind.cov){
							if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=RE.cov[isp,itrans,(Meta$G.transect[isp,itrans]+1):(Meta$M[isp,itrans]/Meta$n.Observers[itrans]),icov]
							else cur.RE=0
							rsamp=switch_sample(n=(Meta$M[isp,itrans]-n.Records[isp,itrans])/Meta$n.Observers[itrans],pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
							Data[isp,itrans,(n.Records[isp,itrans]+1):Meta$M[isp,itrans],Meta$dist.pl+icov]=rep(rsamp,each=Meta$n.Observers[itrans])
						}
					}
          #if(sum(is.na(Data[isp,itrans,,7]))>0)cat(paste("no pop first; iiter ",iiter," isp ",isp," itrans ",itrans))
					
					#update distances, individual covariates for animals that ARE in the population but never observed
					cur.G=Meta$G.transect[isp,itrans]-G.obs[isp,itrans]
					if(cur.G>0){
						#distance
						if(Meta$i.binned==1)dist.star=sample(c(1:Meta$n.bins),cur.G,replace=TRUE,prob=Meta$Bin.length)
						else dist.star=runif(cur.G)
						Cur.dat=Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),]
						if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
						cur.dist=Cur.dat[,Meta$dist.pl]
						X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
						ExpY=X%*%Par$det
						L.old=c(1:length(dist.star))
						Tmp.Y.tilde=Y.tilde[isp,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),itrans]
						for(i in 1:length(dist.star)){
              if(Meta$last.ind)Sigma[offdiag]=Par$cor*(Meta$i.binned*(Meta$n.bins-cur.dist[i])*dist.mult+(1-Meta$i.binned)*cur.dist[i])
							else Sigma[offdiag]=Par$cor*(Meta$i.binned*(cur.dist[i]-1)*dist.mult+(1-Meta$i.binned)*cur.dist[i])
							L.old[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
						}	
						Cur.dat[,Meta$dist.pl]=rep(dist.star,each=Meta$n.Observers[itrans])
						X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
						ExpY=X%*%Par$det
						L.star=L.old
						for(i in 1:length(dist.star)){
						  if(Meta$last.ind)Sigma[offdiag]=Par$cor*(Meta$i.binned*(Meta$n.bins-dist.star[i])*dist.mult+(1-Meta$i.binned)*dist.star[i])
						  else Sigma[offdiag]=Par$cor*(Meta$i.binned*(dist.star[i]-1)*dist.mult+(1-Meta$i.binned)*dist.star[i])
						  L.star[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
						}	
						Acc=(runif(length(L.star))<(L.star/L.old))
						Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl]+rep(Acc,each=Meta$n.Observers[itrans])*rep(dist.star,each=Meta$n.Observers[itrans])
						
						#individual covariates
						if(Meta$n.ind.cov>0){
							for(icov in 1:Meta$n.ind.cov){
								if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=RE.cov[isp,itrans,(G.obs[isp,itrans]+1):Meta$G.transect[isp,itrans],icov]
								else cur.RE=0
								Cov.star=switch_sample(n=cur.G,pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
								Cur.dat=Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),]
								if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
								cur.dist=Cur.dat[,Meta$dist.pl]
								X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
			          #if(iiter==8173 & itrans==21)stop('crap')
                ExpY=X%*%Par$det
								L.old=c(1:length(Cov.star))
								for(i in 1:length(Cov.star)){
								  if(Meta$last.ind)Sigma[offdiag]=Par$cor*(Meta$i.binned*(Meta$n.bins-cur.dist[i])*dist.mult+(1-Meta$i.binned)*cur.dist[i])
								  else Sigma[offdiag]=Par$cor*(Meta$i.binned*(cur.dist[i]-1)*dist.mult+(1-Meta$i.binned)*cur.dist[i])
								  L.old[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
								}	
								Cur.dat[,Meta$dist.pl+icov]=rep(Cov.star,each=Meta$n.Observers[itrans])
								X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
								ExpY=X%*%Par$det
								L.star=L.old
								for(i in 1:length(Cov.star)){
									L.star[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
								}	
								Acc=(runif(length(L.star))<(L.star/L.old))
								Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]+rep(Acc,each=Meta$n.Observers[itrans])*rep(Cov.star,each=Meta$n.Observers[itrans])						
							}
						}
					}
					if(Meta$grps==FALSE)Meta$N.transect[isp,itrans]=Meta$G.transect[isp,itrans]
					else{
					  Meta$N.transect[isp,itrans]=0
					  if(Meta$G.transect[isp,itrans]>0)Meta$N.transect[isp,itrans]=sum(Data[isp,itrans,1:n.Records[isp,itrans],Meta$stacked.names=="Group"])/Meta$n.Observers[itrans]
					}
					#if(sum(is.na(Data[isp,itrans,,7]))>0)cat(paste("not obs first; iiter ",iiter," isp ",isp," itrans ",itrans))
					
					
					#update Y.tilde
					if(n.Records[isp,itrans]>0){
						Cur.dat=Data[isp,itrans,1:n.Records[isp,itrans],]
						if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
						X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
						ExpY=matrix(X%*%Par$det,Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
						Temp.Y.tilde=matrix(Y.tilde[isp,1:n.Records[isp,itrans],itrans],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
						if(Meta$n.Observers[itrans]==1){
							Temp.Y.tilde<-rtruncnorm(Meta$G.transect[isp,itrans],a=ifelse(Cur.dat[,2]==0,-Inf,0),b=ifelse(Cur.dat[,2]==0,0,Inf),ExpY,1)
						}
						else{
							Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
							if(Meta$last.ind)Cor=Par$cor*(Meta$i.binned*(Meta$n.bins-Dist[,1])*dist.mult+(1-Meta$i.binned)*(1-Dist[,1]))
              else Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
							Resp=matrix(Cur.dat[,2],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
							EY1=ExpY[,1]+Cor*(Temp.Y.tilde[,2]-ExpY[,2])
							Temp.Y.tilde[,1] <- rtruncnorm(Meta$G.transect[isp,itrans], a=ifelse(Resp[,1]==0,-Inf,0), b=ifelse(Resp[,1]==0,0,Inf), EY1, sqrt(1-Cor^2))
							EY2=ExpY[,2]+Cor*(Temp.Y.tilde[,1]-ExpY[,1])
							Temp.Y.tilde[,2] <- rtruncnorm(Meta$G.transect[isp,itrans], a=ifelse(Resp[,2]==0,-Inf,0), b=ifelse(Resp[,2]==0,0,Inf), EY2, sqrt(1-Cor^2))
						}
						Y.tilde[isp,1:n.Records[isp,itrans],itrans]=as.vector(t(Temp.Y.tilde))
					}
				}
			}
			if(PROFILE==TRUE){
				cat(paste("Data aug: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
		}

		##### update true species for observed animals ######
		if(Meta$misID){
			#if(iiter%%10==0){ #every 10th iteration 
				#for(isp in 1:Meta$n.species){
					#for(itrans in 1:Meta$n.transects){
						#if(G.obs[isp,itrans]>0){
							#for(iind in 1:G.obs[isp,itrans]){
					
			for(isamp in 1:n.samp.misID){
				isp=sample(Meta$n.species,1,prob=apply(G.obs,1,'sum'))
				itrans=sample(Meta$n.transects,1,prob=G.obs[isp,])
				iind=sample(G.obs[isp,itrans],1)
				Cur.dat=Data[isp,itrans,iind*Meta$n.Observers[itrans]+(1-Meta$n.Observers[itrans]):0,]
				if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
				Sigma=matrix(Par$cor,Meta$n.Observers[itrans],Meta$n.Observers[itrans])
				diag(Sigma)=1
				offdiag=which(Sigma!=1)
				
				if(Meta$detect){
					dist=Cur.dat[1,Meta$dist.pl]
					if(Meta$last.ind)cor=Par$cor*(Meta$i.binned*(Meta$n.bins-dist)*dist.mult+(1-Meta$i.binned)*(1-dist))
					else cor=Par$cor*(Meta$i.binned*(dist-1)*dist.mult+(1-Meta$i.binned)*dist)
          Cur.Y.tilde=Y.tilde[isp,iind*Meta$n.Observers[itrans]+(1-Meta$n.Observers[itrans]):0,itrans]
					X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY=X%*%Par$det
				}
				Other.sp=c(1:Meta$n.species)[-isp]
				if(length(Other.sp)==1)prop.sp=Other.sp
				else prop.sp=sample(Other.sp,1)
				New.dat=Cur.dat
				New.dat[,4]=prop.sp
				if(Meta$detect){
					X=get_mod_matrix(Cur.dat=New.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY.prop=X%*%Par$det
				}
				Obs.species=Cur.dat[,3]
				#abundance component
				mh=log(Lambda.trans[prop.sp,itrans])+log(Meta$G.transect[isp,itrans])-log(Lambda.trans[isp,itrans])-log(Meta$G.transect[prop.sp,itrans]+1)  #verified 4/11/12 & again 6/14/12
				#detection (Y-tilde) component NEED TO PUT BINOMIAL COEFFICIENT IN HERE???
				if(Meta$detect){
					if(Meta$n.Observers[itrans]==1){
						Y.tilde.llik.old=dnorm(Cur.Y.tilde,ExpY,1,log=TRUE)
						Y.tilde.llik.new=dnorm(Cur.Y.tilde,ExpY.prop,1,log=TRUE)
					}
					else{
						Sigma[offdiag]=cor
						Y.tilde.llik.old=dmvnorm(Cur.Y.tilde,ExpY,Sigma,log=TRUE)
						Y.tilde.llik.new=dmvnorm(Cur.Y.tilde,ExpY.prop,Sigma,log=TRUE)						
					}
					mh=mh+Y.tilde.llik.new-Y.tilde.llik.old
				}
				#individual covariate component
				for(icov in 1:Meta$n.ind.cov){
					if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=RE.cov[c(isp,prop.sp),itrans,iind*Meta$n.Observers[itrans],icov]
					else cur.RE=c(0,0)
					cur.cov=Cur.dat[1,Meta$dist.pl+icov]
					mh=mh+switch_pdf(x=cur.cov,pdf=Meta$Cov.prior.pdf[prop.sp,icov],cur.par=Par$Cov.par[prop.sp,,icov],RE=cur.RE)								
					mh=mh-switch_pdf(x=cur.cov,pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,,icov],RE=cur.RE)
				}							
				#misID component		(condition on observed animals)		
				Conf=get_confusion_mat(Cur.dat=Cur.dat,Beta=Par$MisID,misID.mat=Meta$misID.mat,misID.models=Meta$misID.models,misID.symm=Meta$misID.symm,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Levels=Meta$Levels)					
				#Prop.resp=matrix(0,Meta$n.species,n.obs.types)
				for(iobs in 1:Meta$n.Observers[itrans]){
					if(Obs.species[iobs]>0){
						mh=mh-log(Conf[[iobs]][isp,Obs.species[iobs]])+log(Conf[[iobs]][prop.sp,Obs.species[iobs]]) #confusion part
						#Prop.resp[prop.sp,Obs.species[iobs]]=Prop.resp[prop.sp,Obs.species[iobs]]+1
						#Prop.resp[isp,Obs.species[iobs]]=Prop.resp[isp,Obs.species[iobs]]-1
					}
				}
				#Prop.index=rowSums(Prop.resp)
				#for(iindex in 1:Prop.index[prop.sp])mh=mh+log(Multinom.index[itrans,prop.sp]+iindex)
				#for(iobs in 1:n.obs.types)if(Prop.resp[prop.sp,iobs]>0)for(isamp in 1:Prop.resp[prop.sp,iobs])mh=mh-log(Multinom.resp[itrans,prop.sp,iobs]+isamp)
				#for(iindex in 1:(-Prop.index[isp]))mh=mh-log(Multinom.index[itrans,isp]+1-iindex)
				#for(iobs in 1:n.obs.types)if(Prop.resp[isp,iobs]<0)for(isamp in 1:(-Prop.resp[isp,iobs]))mh=mh+log(Multinom.resp[itrans,isp,iobs]+1-isamp)				
				
				if(runif(1)<(exp(mh+log(G.obs[prop.sp,itrans]+1)-log(G.obs[isp,itrans])))){  #accept species change!  (accounting for asymmetry in proposals induced by always proposing 'other' species
					if((n.Records[prop.sp,itrans]+Meta$n.Observers[itrans])<Meta$M[prop.sp,itrans]){						
						#cat(paste("species changed; ind ",iind))
						#1) add data to new target species' data aug matrix
						Data[prop.sp,itrans,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+Meta$n.Observers[itrans]+1):Meta$M[prop.sp,itrans],]=Data[prop.sp,itrans,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1):(Meta$M[prop.sp,itrans]-Meta$n.Observers[itrans]),]	
						Data[prop.sp,itrans,G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1:Meta$n.Observers[itrans],]=New.dat
						#2) remove entry from old species; shift data aug array down
						Data[isp,itrans,((iind-1)*(Meta$n.Observers[itrans])+1):(Meta$M[isp,itrans]-Meta$n.Observers[itrans]),]=Data[isp,itrans,((iind-1)*(Meta$n.Observers[itrans])+1+Meta$n.Observers[itrans]):Meta$M[isp,itrans],]
						#3) add Y.tilde to new target species Y.tilde matrix
						if(Meta$detect){
							Y.tilde[prop.sp,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+Meta$n.Observers[itrans]+1):Meta$M[prop.sp,itrans],itrans]=Y.tilde[prop.sp,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1):(Meta$M[prop.sp,itrans]-Meta$n.Observers[itrans]),itrans]	
							Y.tilde[prop.sp,G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1:Meta$n.Observers[itrans],itrans]=Cur.Y.tilde	
							#4) remove Y.tilde from current species' Y.tilde matrix
							Y.tilde[isp,((iind-1)*(Meta$n.Observers[itrans])+1):(Meta$M[isp,itrans]-Meta$n.Observers[itrans]),itrans]=Y.tilde[isp,((iind-1)*(Meta$n.Observers[itrans])+1+Meta$n.Observers[itrans]):Meta$M[isp,itrans],itrans]
						}
						#5) adjust dimensions of things
						G.obs[isp,itrans]=G.obs[isp,itrans]-1
						Meta$G.transect[isp,itrans]=Meta$G.transect[isp,itrans]-1
						n.Records[isp,itrans]=n.Records[isp,itrans]-Meta$n.Observers[itrans]
						G.obs[prop.sp,itrans]=G.obs[prop.sp,itrans]+1
						Meta$G.transect[prop.sp,itrans]=Meta$G.transect[prop.sp,itrans]+1
						n.Records[prop.sp,itrans]=n.Records[prop.sp,itrans]+Meta$n.Observers[itrans]
						if(Meta$grps==FALSE){
							Meta$N.transect[isp,itrans]=Meta$G.transect[isp,itrans]
							Meta$N.transect[prop.sp,itrans]=Meta$G.transect[prop.sp,itrans]
						}
						else{
							Meta$N.transect[isp,itrans]=Meta$N.transect[isp,itrans]-Cur.dat[1,grp.pl]
							Meta$N.transect[prop.sp,itrans]=Meta$N.transect[prop.sp,itrans]+Cur.dat[1,grp.pl]
						}
						#6) update multinomial index computation stuff
						#Multinom.index[itrans,]=Multinom.index[itrans,]+Prop.index
						#Multinom.resp[itrans,,]=Multinom.resp[itrans,,]+Prop.resp
					}
				}		
			}#}}}}
			if(PROFILE==TRUE){
				cat(paste("True species: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
		
			
			##### update misID parameters 

				#note: can no longer do by species since misID.symm=TRUE means some parameters affect multiple species         
			  Cur.dat=stack_data_misID(Data=Data,G.obs=G.obs,g.tot.obs=g.tot.obs,n.Observers=Meta$n.Observers,n.transects=Meta$n.transects,n.species=Meta$n.species,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind)
        if(sum(Cur.dat[,3]==0)>0)Cur.dat=Cur.dat[-which(Cur.dat[,3]==0),]  #eliminate records where animals were missed
         
				Conf=get_confusion_mat(Cur.dat=Cur.dat,Beta=Par$MisID,misID.mat=Meta$misID.mat,misID.models=Meta$misID.models,misID.symm=Meta$misID.symm,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Levels=Meta$Levels)					
				Probs=rep(0,nrow(Cur.dat))
				for(itmp in 1:nrow(Cur.dat))Probs[itmp]=Conf[[itmp]][Cur.dat[itmp,4],Cur.dat[itmp,3]]
				
				logL.old=sum(log(Probs))  #categorical distribution
				
				Temp.par=Par$MisID
				for(ipl in 1:max(Meta$misID.mat)){  #now, update parameter values
					for(ipar in 1:Meta$N.par.misID[ipl]){
						beta.old=Par$MisID[[ipl]][ipar]
						beta.star=beta.old+rnorm(1,0,Control$MH.misID[ipl,ipar])
						Temp.par[[ipl]][ipar]=beta.star
						Conf.new=get_confusion_mat(Cur.dat=Cur.dat,Beta=Temp.par,misID.mat=Meta$misID.mat,misID.models=Meta$misID.models,misID.symm=Meta$misID.symm,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Levels=Meta$Levels)					
						#make sure Conf.new has [isp,isp] entries larger than [isp,other sp] entries to prevent label swapping	
						Conf.sums=0*Conf[[1]]
						flag=0
						for(isp in 1:Meta$n.species){
							for(icol in 1:n.obs.types)Conf.sums[isp,icol]=sum(sapply(Conf.new,function(x,isp,icol)x[isp,icol],isp=isp,icol=icol))
							if(Conf.sums[isp,isp]!=max(Conf.sums[isp,]))flag=1
						}
						if(flag==0){ 
							for(itmp in 1:nrow(Cur.dat))Probs[itmp]=Conf.new[[itmp]][Cur.dat[itmp,4],Cur.dat[itmp,3]]
							logL.new=sum(log(Probs))  #categorical distribution
							if(runif(1)<exp(logL.new-logL.old+dnorm(beta.star,Prior.pars$misID.mu[[ipl]][ipar],Prior.pars$misID.sd[[ipl]][ipar],log=1)-dnorm(beta.old,Prior.pars$misID.mu[[ipl]][ipar],Prior.pars$misID.sd[[ipl]][ipar],log=1))){
								Par$MisID[[ipl]][ipar]=beta.star
								Accept$MisID[[ipl]][ipar]=Accept$MisID[[ipl]][ipar]+1
								logL.old=logL.new
							}
							else Temp.par[[ipl]][ipar]=beta.old
						}
						else Temp.par[[ipl]][ipar]=beta.old
					}
				}
			}
			
			if(PROFILE==TRUE){
				cat(paste("misID pars: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}		
		
		
		if(Meta$detect){
			###############       update detection process parameters       ##############
			# First, assemble stacked adjusted Response, X matrices across all transects; 
			#basic form of response is Ytilde[obs1]-cor*Ytilde[obs2]
			#adjusted DM rows are of form X[obs1]-cor*X[obs2]
			X.beta=matrix(NA,1,n.beta.det)
			Y.beta=NA
			Cor=NA
			for(isp in 1:Meta$n.species){
				GT0=which(Meta$G.transect[isp,]>0)
				n.gt0=length(GT0)
				if(n.gt0>0){
					for(itrans in 1:n.gt0){
						Cur.dat=Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],]
						if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
						Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1]
						X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,Meta$n.Observers[GT0[itrans]],Meta$G.transect[isp,GT0[itrans]]))
						if(Meta$n.Observers[GT0[itrans]]==2){
              if(Meta$last.ind)Tmp.cor=Par$cor*(Meta$i.binned*(Meta$n.bins-Dist)*dist.mult+(1-Meta$i.binned)*(1-Dist))
						  else Tmp.cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
						  Cor=c(Cor,rep(Tmp.cor,2))  #assemble vector of correlation parameters for each observation
              X.beta=rbind(X.beta,t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
						  Y.temp=matrix(Y.tilde[isp,1:n.Records[isp,GT0[itrans]],GT0[itrans]],Meta$G.transect[isp,GT0[itrans]],2,byrow=TRUE)
						  Y.beta=c(Y.beta,Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
						}
						else{
							X.beta=rbind(X.beta,t(X.temp[,1,]))
							Y.beta=c(Y.beta,Y.tilde[isp,1:n.Records[isp,GT0[itrans]],GT0[itrans]])
							Cor=c(Cor,rep(0,n.Records[isp,GT0[itrans]]))
						}
					}
				}
			}
			X.beta=X.beta[-1,]
			Y.beta=Y.beta[-1]
			Cor=Cor[-1]
			#now use basic matrix equations from Gelman '04 book (14.11 14.12) to update beta parms
			Sig.inv=1/(1-Cor^2)  #for use in eq. 14.11, 14.12 of Gelman et al.; don't need a matrix since it would be diagonal
			V.inv <- crossprod(X.beta,Sig.inv*X.beta) 	
			M.z <- solve(V.inv, crossprod(X.beta,Sig.inv*Y.beta))
			Par$det <- M.z + solve(chol(V.inv), rnorm(n.beta.det,0,1))
			
			if(PROFILE==TRUE){
				cat(paste("Detection pars: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			
			
			#update correlation parameter for detection process (if applicable)
			if(Meta$point.ind==1){
				cor.star=Par$cor+runif(1,-Control$MH.cor,Control$MH.cor)
				if(cor.star>max(-0.95,-0.95*(1-(Meta$last.ind==FALSE & Meta$cor.const==TRUE))) & cor.star<min(0.95,0.95*(1-(Meta$last.ind==TRUE & Meta$cor.const==TRUE)))){
					Delta1=rep(NA,sum(Meta$G.transect[,which(Meta$n.Observers==2)]))
					Delta2=Delta1
					Dist=Delta1
					counter=1
					for(isp in 1:Meta$n.species){
						I.gt.one=which(Meta$n.Observers>1 & Meta$G.transect[isp,]>0)
						n.gt.one=length(I.gt.one)
						if(n.gt.one>0){
							for(itrans in 1:n.gt.one){
								Cur.dat=Data[isp,I.gt.one[itrans],1:n.Records[isp,I.gt.one[itrans]],]
								if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
								Dist[counter:(counter+Meta$G.transect[isp,I.gt.one[itrans]]-1)]=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[isp,I.gt.one[itrans]],Meta$n.Observers[I.gt.one[itrans]],byrow=TRUE)[,1]
								X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[isp,I.gt.one[itrans]]))
								Y.temp=matrix(Y.tilde[isp,1:n.Records[isp,I.gt.one[itrans]],I.gt.one[itrans]],Meta$G.transect[isp,I.gt.one[itrans]],2,byrow=TRUE)			
								Delta1[counter:(counter+Meta$G.transect[isp,I.gt.one[itrans]]-1)]=Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det)
								Delta2[counter:(counter+Meta$G.transect[isp,I.gt.one[itrans]]-1)]=Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det)
								counter=counter+Meta$G.transect[isp,I.gt.one[itrans]]
							}
						}
					}
          if(Meta$last.ind)Cor=Par$cor*(Meta$i.binned*(Meta$n.bins-Dist)*dist.mult+(1-Meta$i.binned)*(1-Dist))
					else Cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
					logP.old=-.5*(sum(log(1-Cor^2))+sum((Delta1^2+Delta2^2-2*Cor*Delta1*Delta2)/(1-Cor^2)))
					if(Meta$last.ind)Cor=cor.star*(Meta$i.binned*(Meta$n.bins-Dist)*dist.mult+(1-Meta$i.binned)*(1-Dist))
					else Cor=cor.star*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
          logP.new=-.5*(sum(log(1-Cor^2))+sum((Delta1^2+Delta2^2-2*Cor*Delta1*Delta2)/(1-Cor^2)))
					if(runif(1)<exp(logP.new-logP.old)){
						Par$cor=cor.star
						Accept$cor=Accept$cor+1
					}				
				}
			}
			
			if(PROFILE==TRUE){
				cat(paste("Correlation: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
		}
		#update parameters of individual covariate distributions (if fixed=0)
		for(isp in 1:Meta$n.species){
      if(sum(Meta$G.transect[isp,])>0){
        GT0=which(Meta$G.transect[isp,]>0)
        n.gt0=length(GT0)
        for(icov in 1:Meta$n.ind.cov){
          if(Meta$Cov.prior.fixed[isp,icov]==0){
            if(Meta$Cov.prior.pdf[isp,icov]=="normal")cat("\n Warning: hyper-priors not yet implemented for normal dist. \n")
            if(Meta$Cov.prior.pdf[isp,icov]=="poisson"){
              Cur.cov=matrix(Data[isp,GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
              if(n.gt0>1){
                for(itrans in 2:n.gt0){
                  Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
                }
              }
              Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="pois1"){
              Cur.cov=matrix(Data[isp,GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
              if(n.gt0>1){
                for(itrans in 2:n.gt0){
                  Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
                }
              }
              Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)-length(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="poisson_ln" | Meta$Cov.prior.pdf[isp,icov]=="pois1_ln"){
              Cur.cov=matrix(Data[isp,GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
              Cur.RE=RE.cov[isp,GT0[1],1:Meta$G.transect[isp,GT0[1]],icov]
              if(n.gt0>1){
                for(itrans in 2:n.gt0){
                  Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
                  Cur.RE=c(Cur.RE,RE.cov[isp,GT0[itrans],1:Meta$G.transect[isp,GT0[itrans]],icov])
                }
              }
              Cur.cov=Cur.cov-(Meta$Cov.prior.pdf[isp,icov]=="pois1_ln")
              #1) update theta
              par.star=Par$Cov.par[isp,1,icov]+runif(1,-0.05,0.05)
              sum.y=sum(Cur.cov)
              sum.yZ=sum(Cur.cov*Cur.RE)
              log.post.new=par.star*sum.y-sum(exp(par.star+Par$Cov.par[isp,2,icov]*Cur.RE))+dnorm(par.star,Meta$Cov.prior.parms[isp,1,icov],Meta$Cov.prior.parms[isp,2,icov],log=1)
              log.post.old=Par$Cov.par[isp,1,icov]*sum.y-sum(exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Cur.RE))+dnorm(Par$Cov.par[isp,1,icov],Meta$Cov.prior.parms[isp,1,icov],Meta$Cov.prior.parms[isp,2,icov],log=1)
              if(runif(1)<exp(log.post.new-log.post.old))Par$Cov.par[isp,1,icov]=par.star
              #2) update sigma
              par.star=Par$Cov.par[isp,2,icov]+runif(1,-.01,.01)
              if(par.star>0 & par.star<Meta$Cov.prior.parms[isp,3,icov]){
                log.post.new=par.star*sum.yZ-sum(exp(Par$Cov.par[isp,1,icov]+par.star*Cur.RE))
                log.post.old=Par$Cov.par[isp,2,icov]*sum.yZ-sum(exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Cur.RE))
                if(runif(1)<exp(log.post.new-log.post.old))Par$Cov.par[isp,2,icov]=par.star
              }
              #3) update random effects		
              for(itrans in 1:n.gt0){
                #animals currently in population
                Cur.cov=matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1]-(Meta$Cov.prior.pdf[isp,icov]=="pois1_ln")						
                Cur.RE=RE.cov[isp,GT0[itrans],1:Meta$G.transect[isp,GT0[itrans]],icov]
                Prop=Cur.RE+runif(length(Cur.RE),-.1,.1)
                LogPost.new=Cur.cov*Par$Cov.par[isp,2,icov]*Prop-exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Prop)+dnorm(Prop,0,1,log=1)
                LogPost.old=Cur.cov*Par$Cov.par[isp,2,icov]*Cur.RE-exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Cur.RE)+dnorm(Cur.RE,0,1,log=1)
                Acc=(runif(length(Cur.RE))<(exp(LogPost.new-LogPost.old)))
                RE.cov[isp,GT0[itrans],1:Meta$G.transect[isp,GT0[itrans]],icov]=Acc*Prop+(1-Acc)*Cur.RE
              }
              #animals currently not in population
              for(itrans in 1:Meta$n.transects){
                RE.cov[isp,itrans,(Meta$G.transect[isp,itrans]+1):Meta$M[isp,itrans],icov]=rnorm(Meta$M[isp,itrans]-Meta$G.transect[isp,itrans],0,1)
              }
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="multinom"){
              Cur.cov=matrix(Data[isp,1:GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
              if(n.gt0>1){
                for(itrans in 2:n.gt0){
                  Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
                }
              }
              Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov]=rdirichlet(1,Meta$Cov.prior.parms[isp,1:Meta$Cov.prior.n[isp,icov],icov]+tabulate(factor(Cur.cov)))
            }
          }
        }
      }
		}
    #if(is.na(Par$Cov.par[1,1,1])){
    #  cat("par first, iter ")
    #  cat(iiter)
    #}
		if(PROFILE==TRUE){
			cat(paste("Ind cov pars: ", (Sys.time()-st),'\n'))
			st=Sys.time()
		}
		
		if(adapt==TRUE){
			if(iiter%%100==0){
				if(Accept$cor<30)Control$MH.cor=Control$MH.cor*.95
				if(Accept$cor>40)Control$MH.cor=Control$MH.cor*1.053
				if(Meta$misID){
					for(ipl in 1:length(Meta$N.par.misID)){
						for(ipar in 1:Meta$N.par.misID[ipl]){
							if(Accept$MisID[[ipl]][ipar]<30)Control$MH.misID[ipl,ipar]=Control$MH.misID[ipl,ipar]*0.95
							if(Accept$MisID[[ipl]][ipar]>40)Control$MH.misID[ipl,ipar]=Control$MH.misID[ipl,ipar]*1.053
						}
					}
				}
				for(ipar in 1:Meta$n.species){
					for(i in 1:n.unique)
					if(Accept$Nu[ipar,i]<30)Control$MH.nu[ipar,i]=Control$MH.nu[ipar,i]*.95
					if(Accept$Nu[ipar,i]>40)Control$MH.nu[ipar,i]=Control$MH.nu[ipar,i]*1.053
				}
				for(ipar in 1:length(Accept$N)){
					if(Accept$N[ipar]<30)Control$RJ.N[ipar]=max(1,Control$RJ.N[ipar]-1)
					if(Accept$N[ipar]>40)Control$RJ.N[ipar]=Control$RJ.N[ipar]+1
				}
				Accept$cor=0
				Accept$Hab=Accept$Hab*0
				Accept$Nu=Accept$Nu*0
				Accept$N=Accept$N*0
				Accept$MisID=lapply(Accept$MisID,function(x)x*0)
			}
		}
		
		
		#store results if applicable
		if(iiter>Control$burnin & iiter%%Control$thin==0){
			MCMC$cor[(iiter-Control$burnin)/Control$thin]=Par$cor
			MCMC$Det[(iiter-Control$burnin)/Control$thin,]=Par$det
			if(Meta$misID)for(i in 1:length(Meta$N.par.misID))MCMC$MisID[[i]][,(iiter-Control$burnin)/Control$thin]=Par$MisID[[i]]
			for(isp in 1:Meta$n.species){
				MCMC$G[isp,(iiter-Control$burnin)/Control$thin,]=Par$G[isp,]
				MCMC$N[isp,(iiter-Control$burnin)/Control$thin,]=Par$N[isp,]
				MCMC$N.tot[isp,(iiter-Control$burnin)/Control$thin]=sum(Par$N[isp,])
				MCMC$Hab.pois[isp,(iiter-Control$burnin)/Control$thin,]=Par$hab.pois[isp,]
				MCMC$tau.eta.pois[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eta.pois[isp]
        if(Meta$ZIP){
				  MCMC$Hab.bern[isp,(iiter-Control$burnin)/Control$thin,]=Par$hab.bern[isp,]
				  MCMC$tau.eta.bern[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eta.bern[isp]
        }
				MCMC$tau.nu[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.nu[isp]
				MCMC$Cov.par[isp,(iiter-Control$burnin)/Control$thin,]=Par$Cov.par[isp,,]
				Obs.N[isp,(iiter-Control$burnin)/Control$thin,]=Meta$N.transect[isp,]
				Temp.G=Meta$Area.hab[Meta$Mapping]*Meta$Area.trans*exp(rnorm(Meta$n.transects,(DM.hab.pois[[isp]]%*%Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]+Par$Eta.pois[isp,])[Meta$Mapping],sqrt(1/Par$tau.nu[isp])))
        if(Meta$ZIP)Temp.G=Temp.G*rbern(Meta$n.transects,pnorm((DM.hab.bern[[isp]]%*%Par$hab.bern[isp,1:Meta$N.hab.bern.par[isp]]+Par$Eta.bern[isp,])[Meta$Mapping]))
        Pred.N[isp,(iiter-Control$burnin)/Control$thin,]=Temp.G+rpois(Meta$n.transects,grp.lam[isp]*Temp.G)	
			}
       #posterior predictions of detection data given nu, detection & misclassification parameters
			if(Meta$post.loss){ #calculate observed counts of different detection types, initialize prediction arrays
        Sigma=diag(2)
        Cur.lambda=exp(Par$Nu)*Meta$Area.hab[Meta$Mapping]
        for(isp in 1:Meta$n.species)Cur.lambda[isp,]=Cur.lambda[isp,]*Meta$Area.trans
        Cur.G=matrix(rpois(Meta$n.species*Meta$n.transects,Cur.lambda),Meta$n.species,Meta$n.transects)		
        if(Meta$ZIP)Cur.G=Z*Cur.G         
        for(itrans in 1:Meta$n.transects){
          for(isp in 1:Meta$n.species){
            if(Cur.G[isp,itrans]>0){
              Cur.dat=matrix(0,Cur.G[isp,itrans]*Meta$n.Observers[itrans],dim(Data)[4])
              Cur.dat[,3]=isp
              #fill observer
              if(Meta$n.Observers[itrans]==1)Cur.dat[,1]=Data[1,itrans,1,1]
              else Cur.dat[,1]=Data[1,itrans,1:2,1]
              #fill observer covariates
              if(n.obs.cov>0){
                for(ipl in 4:(3+n.obs.cov)){
                  Cur.dat[,ipl]=rep(Data[isp,itrans,1:Meta$n.Observers[itrans],ipl],Cur.G[isp,itrans])
                }
              }
              #sample distance
              if(Meta$i.binned==1)Cur.dat[,Meta$dist.pl]=rep(sample(c(1:Meta$n.bins),size=Cur.G[isp,itrans],prob=Meta$Bin.length,replace=TRUE),each=Meta$n.Observers[itrans])
              else Cur.dat[,Meta$dist.pl]=rep(runif(Cur.G[isp,itrans]),each=Meta$n.Observers[itrans])
              #sample from individual covariate distributions 
              if(Meta$n.ind.cov>0){
                for(icov in 1:Meta$n.ind.cov){
                  if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=rep(rnorm(Cur.G[isp,itrans],0,1),each=Meta$n.Observers[itrans])
                  else cur.RE=0
                  rsamp=switch_sample(n=Cur.G[isp,itrans],pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
                  Cur.dat[,Meta$dist.pl+icov]=rep(rsamp,each=Meta$n.Observers[itrans])
                }
              }
              if(Meta$detect)X.temp=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
              if(Meta$n.Observers[itrans]==1){ #in this case, univariate detection; just fill first column of Pred.det
                if(Meta$detect==FALSE)Cur.dat[,2]=1
                else Cur.dat[,2]=(rnorm(Cur.G[isp,itrans],X.temp%*%Par$det,1)>0) #probit detection model
                if(sum(Cur.dat[,2])>0){ #misID model (if applicable)
                  if(Meta$misID){
                    Det.ind=which(Cur.dat[,2]==1)
                    Conf=get_confusion_mat(Cur.dat=matrix(Cur.dat[Det.ind,],nrow=length(Det.ind)),Beta=Par$MisID,misID.mat=Meta$misID.mat,misID.models=Meta$misID.models,misID.symm=Meta$misID.symm,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Levels=Meta$Levels)  				
                    for(iind in 1:length(Det.ind)){
                      cur.sp=sample(1:ncol(Conf[[iind]]),1,prob=Conf[[iind]][isp,])
                      Pred.det[(iiter-Control$burnin)/Control$thin,itrans,cur.sp+1,1]=Pred.det[(iiter-Control$burnin)/Control$thin,itrans,cur.sp+1,1]+1                    
                    }
                  }
                  else{
                    Det.ind=which(Cur.dat[,2]==1)
                    for(iind in 1:length(Det.ind)){
                      cur.sp=Cur.dat[Det.ind[iind],3]               
                      Pred.det[(iiter-Control$burnin)/Control$thin,itrans,cur.sp+1,1]=Pred.det[(iiter-Control$burnin)/Control$thin,itrans,cur.sp+1,1]+1                                        
                    }
                  }
                }
              }              
              else{  #in this case, bivariate detection
                if(Meta$detect)XB=X.temp%*%Par$det
                for(iind in 1:Cur.G[isp,itrans]){
                  if(Meta$point.ind & Meta$detect){
                    cur.dist=Cur.dat[iind*2,Meta$dist.pl]
                    if(Meta$last.ind)Sigma[offdiag]=Par$cor*(Meta$i.binned*(Meta$n.bins-cur.dist)*dist.mult+(1-Meta$i.binned)*cur.dist)
                    else Sigma[offdiag]=Par$cor*(Meta$i.binned*(cur.dist-1)*dist.mult+(1-Meta$i.binned)*cur.dist)
                  }    
                  if(Meta$detect)Cur.det=(rmvnorm(1,XB[(iind*2-1):(iind*2)],Sigma)>0)  #bivariate normal detection
                  else Cur.det=c(1,1)
                  if(sum(Cur.det)>0){
                    if(Meta$misID){
                      Cur.obs=rep(0,2)
                      for(iobs in 1:2){
                        if(Cur.det[iobs]==1){ #only model misID for detections
                          Conf=get_confusion_mat(Cur.dat=matrix(Cur.dat[iind*2-2+iobs,],nrow=1),Beta=Par$MisID,misID.mat=Meta$misID.mat,misID.models=Meta$misID.models,misID.symm=Meta$misID.symm,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Levels=Meta$Levels)    			
                          Cur.obs[iobs]=sample(1:ncol(Conf[[1]]),1,prob=Conf[[1]][isp,])
                        }
                      }
                    }
                    else Cur.obs=rep(Cur.dat[iind*2,3],2)  
                    Pred.det[(iiter-Control$burnin)/Control$thin,itrans,Cur.obs[1]+1,Cur.obs[2]+1]=Pred.det[(iiter-Control$burnin)/Control$thin,itrans,Cur.obs[1]+1,Cur.obs[2]+1]+1                    
                  }
                }
              } 
            }
          }
			  }
			}			
		}
		
		if(iiter==100){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/100
			ttc <- round((cur.iter-100)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}	
		if((iiter%%1000)==1)cat(paste('iteration ', iiter,' of ',cur.iter,' completed \n'))
	}
	cat(paste('\n total elapsed time: ',difftime(Sys.time(),st,units="mins"),' minutes \n'))

	Post=list(N=MCMC$N,G=MCMC$G)
	#convert Out$MCMC into mcmc object for use with coda, S3 methods
	Hab.pois.names=vector("list",Meta$n.species)
	for(isp in 1:Meta$n.species)Hab.pois.names[[isp]]=colnames(DM.hab.pois[[isp]])
  if(Meta$ZIP){
    Hab.bern.names=vector("list",Meta$n.species)
    for(isp in 1:Meta$n.species)Hab.bern.names[[isp]]=colnames(DM.hab.bern[[isp]])
  }
  Det.names=colnames(get_mod_matrix(Cur.dat=matrix(Data[1,1,1,],nrow=1),Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels))
	Cov.names=vector("list",Meta$n.ind.cov)
	Cov.par.n=0
	if(Meta$n.ind.cov>0){
		for(icov in 1:Meta$n.ind.cov){
			Par.name=switch(Meta$Cov.prior.pdf[icov],pois1_ln=c("mean.minus.1","sd"),poisson_ln=c("mean","sd"),multinom=paste("prop.cell.",c(1:(Meta$Cov.prior.n[isp,icov]-1)),sep=''),normal="mean",pois1="mean.minus.1",poisson="mean")
			Cov.names[[icov]]=paste(Meta$stacked.names[Meta$dist.pl+icov],".",Par.name,sep='')
			Cov.par.n=Cov.par.n+length(Cov.names[[icov]])
		}
	}
	MisID.names=NULL
	if(Meta$misID==TRUE){
		MisID.names=vector("list",max(Meta$misID.mat))
		for(imod in 1:max(Meta$misID.mat)){
			MisID.names[[imod]]=colnames(get_mod_matrix(Cur.dat=matrix(Data[1,1,1,],nrow=1),stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[imod]],Levels=Meta$Levels))		
		}
	}	
  if(1==1){
    if(Meta$ZIP){
      N.hab.bern.par=Meta$N.hab.bern.par
      Hab.bern.names=Hab.bern.names
    }
    else{
      N.hab.bern.par=NA
      Hab.bern.names=NA
    }
  }
	MCMC=convert.HDS.to.mcmc(MCMC=MCMC,N.hab.pois.par=Meta$N.hab.pois.par,N.hab.bern.par=N.hab.bern.par,Cov.par.n=Cov.par.n,Hab.pois.names=Hab.pois.names,Hab.bern.names=Hab.bern.names,Det.names=Det.names,Cov.names=Cov.names,MisID.names=MisID.names,N.par.misID=Meta$N.par.misID,misID.mat=Meta$misID.mat,fix.tau.nu=Meta$fix.tau.nu,misID=Meta$misID,spat.ind=Meta$spat.ind,point.ind=Meta$point.ind)
	Out=list(Post=Post,MCMC=MCMC,Accept=Accept,Control=Control,Obs.N=Obs.N,Pred.N=Pred.N,Obs.det=Obs.det,Pred.det=Pred.det)
	Out
}

