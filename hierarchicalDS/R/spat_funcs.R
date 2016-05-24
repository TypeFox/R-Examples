if(getRversion() >= "2.15.1")  utils::globalVariables(c("Easting", "Northing"))

#' function to sample from a specified probability density function
#' @param n number of samples desired
#' @param pdf probability density function (pois1, poisson, normal, unif.disc, unif.cont)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @param RE random effects, if present
#' @return a vector of length n samples from the desired distribution 
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_sample<-function(n,pdf,cur.par,RE){
  switch(pdf,
         pois1=rpois(n,cur.par[1])+1,
         poisson=rpois(n,cur.par[1]),
         pois1_ln=rpois(n,exp(cur.par[1]+cur.par[2]*RE))+1,
         poisson_ln=rpois(n,exp(cur.par[1]+cur.par[2]*RE)),
         normal=rnorm(n,cur.par[1],cur.par[2]),
         unif.disc=sample(cur.par[1]:cur.par[2],n,replace=TRUE),
         unif.cont=runif(n,cur.par[1],cur.par[2]),
         multinom=sample(c(1:length(cur.par)),n,replace=TRUE,prob=cur.par)
  )
}

#' function to sample from hyperpriors of a specified probability density function; note that
#' initial values for sigma of lognormal random effects are fixed to a small value (0.05) to
#' prevent numerical errors
#' @param pdf probability density function (pois1, poisson, normal, unif.disc, unif.cont)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @return a vector of length n samples from the desired distribution 
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_sample_prior<-function(pdf,cur.par){
  #require(mc2d)
  switch(pdf,
         pois1=rgamma(1,cur.par[1],cur.par[2]),
         poisson=rgamma(1,cur.par[1],cur.par[2]),
         pois1_ln=c(rnorm(1,cur.par[1],cur.par[2]),0.05),
         poisson_ln=c(rnorm(1,cur.par[1],cur.par[2]),0.05),
         multinom=rdirichlet(1,cur.par)
  )
}

#' function to calculate the joint pdf for a sample of values from one of a number of pdfs
#' @param x values to be evaluated
#' @param pdf probability density function (pois1, poisson, pois1_ln, poisson_ln, normal, multinom)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @param RE random effects, if present
#' @return total log likelihood of points
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_pdf<-function(x,pdf,cur.par,RE){
  switch(pdf,
         pois1=sum(dpois(x-1,cur.par[1],log=1)),
         poisson=sum(dpois(x,cur.par[1],log=1)),
         pois1_ln=sum(dpois(x-1,exp(cur.par[1]+cur.par[2]*RE),log=1)),
         poisson_ln=sum(dpois(x,exp(cur.par[1]+cur.par[2]*RE),log=1)),
         normal=sum(dnorm(x,cur.par[1],cur.par[2],log=1)),
         multinom=sum(log(cur.par[x]))
  )
}

#' function to stack data (going from three dimensional array to a two dimensional array including only "existing" animals
#' @param Data three-d dataset
#' @param Obs.transect	current number of observations of animals in each transect (vector)
#' @param n.transects	number of transects
#' @param stacked.names	column names for new stacked dataset
#' @param factor.ind	a vector of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
#' @return a stacked dataset
#' @export
#' @keywords stack data
#' @author Paul B. Conn
stack_data<-function(Data,Obs.transect,n.transects,stacked.names,factor.ind){
  #convert from "sparse" 3-d data augmentation array to a rich 2-d dataframe for updating beta parameters 
  if(n.transects==1)Stacked=Data
  else{
    Stacked=as.data.frame(Data[1,1:2,])
    for(itrans in 1:n.transects){
      if(Obs.transect[itrans]>0)Stacked=rbind(Stacked,Data[itrans,1:Obs.transect[itrans],])
    }
    Stacked=Stacked[-c(1,2),]
  }
  colnames(Stacked)=stacked.names	#gotta reestablish variable type since 3-d array doesn't hold it
  factor.cols=which(factor.ind[stacked.names]==TRUE) 
  if(length(factor.cols)>0){
    for(icol in 1:length(factor.cols)){
      Stacked[,factor.cols[icol]]=as.factor(Stacked[,factor.cols[icol]])
    }
  }
  Stacked
}

#' function to stack data for midID updates (going from four dimensional array to a two dimensional array including observed groups
#' @param Data 4-d dataset
#' @param G.obs matrix giving the total numer of groups observed at least once by species and transect
#' @param g.tot.obs  total number of observations for animals seen at least once
#' @param n.Observers vector giving number of observers per transect
#' @param n.transects  number of transects
#' @param n.species number of species
#' @param stacked.names  column names for new stacked dataset
#' @param factor.ind  a vector of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
#' @return a stacked dataset (in matrix form)
#' @export
#' @keywords stack data
#' @author Paul B. Conn
stack_data_misID<-function(Data,G.obs,g.tot.obs,n.Observers,n.transects,n.species,stacked.names,factor.ind){
  #convert from "sparse" 4-d data augmentation array to a rich 2-d dataframe for updating misID parameters 
  if(n.transects==1 & n.species==1)Stacked=Data[1,1,,]
  else{
    G.tot.obs=G.obs
    Stacked=matrix(0,g.tot.obs,length(Data[1,1,1,]))
    ipl=1
    for(isp in 1:n.species){
      G.tot.obs[isp,]=G.obs[isp,]*n.Observers
      for(itrans in 1:n.transects){
        if(G.obs[isp,itrans]>0)Stacked[ipl:(ipl+G.tot.obs[isp,itrans]-1),]=Data[isp,itrans,1:G.tot.obs[isp,itrans],]
        ipl=ipl+G.tot.obs[isp,itrans]
      }
    }
  }
  Stacked
}


#' function to produce a design matrix given a dataset and user-specified formula object
#' @param Cur.dat 	current dataset
#' @param stacked.names	column names for current dataset
#' @param factor.ind	a list of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
#' @param Det.formula	a formula object
#' @param Levels	A list object giving the number of levels for factor variables
#' @return a design matrix
#' @export
#' @keywords model matrix
#' @author Paul B. Conn
get_mod_matrix<-function(Cur.dat,stacked.names,factor.ind,Det.formula,Levels){
  Cur.dat=as.data.frame(Cur.dat)
  colnames(Cur.dat)=stacked.names
  factor.cols=which(factor.ind[stacked.names]==TRUE)
  if(length(factor.cols)>0){
    for(icol in 1:length(factor.cols)){
      Cur.dat[,factor.cols[icol]]=eval(parse(text=paste('factor(Cur.dat[,factor.cols[icol]],levels=Levels$',names(factor.cols)[icol],')',sep='')))
    }
  }
  DM=model.matrix(Det.formula,data=Cur.dat)
  DM
}

#' generate initial values for MCMC chain if not already specified by user
#' @param DM.hab 	design matrix for habitat model
#' @param DM.det	design matrix for detection model
#' @param G.transect a vector of the number of groups of animals in area covered by each transect		
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param point.ind	is point independence assumed (TRUE/FALSE)
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  pois1 parameter for group size
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits<-function(DM.hab,DM.det,G.transect,Area.trans,Area.hab,Mapping,point.ind,spat.ind,grp.mean){		
  Par=list(det=rnorm(ncol(DM.det),0,1),hab=rep(0,ncol(DM.hab)),cor=ifelse(point.ind,runif(1,0,.8),0),
           Nu=log(max(G.transect)/mean(Area.trans)*exp(rnorm(length(Area.hab)))),Eta=rnorm(length(Area.hab)),
           tau.eta=runif(1,0.5,2),tau.nu=runif(1,0.5,2))
  Par$hab[1]=mean(G.transect)/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(1,0,1))
  Par$G=round(exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu))))
  Par$N=Par$G+rpois(length(Par$G),grp.mean*Par$G)
  if(spat.ind==1)Par$Eta=0*Par$Eta
  Par
}

#' generate initial values for misID model if not already specified by user
#' @param DM.hab.pois 	a list of design matrices for the Poisson habitat model (elements are named sp1,sp2, etc.)
#' @param DM.hab.bern   If a hurdle model, a list of design matrices for the Bernoulli habitat model (elements are named sp1,sp2, etc.) (NULL if not hurdle)
#' @param DM.det	design matrix for detection model
#' @param N.hab.pois.par  vector giving number of parameters in the Poisson habitat model for each species
#' @param N.hab.bern.par  vector giving number of parameters in the Bernoulli habitat model for each species (NULL if not hurdle)
#' @param G.transect a matrix of the number of groups of animals in area covered by each transect; each row gives a separate species		
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param point.ind	is point independence assumed (TRUE/FALSE)
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  a vector giving the pois1 parameter for group size (one entry for each species)
#' @param misID    if TRUE, indicates that misidentification is incorporated into modeling
#' @param misID.mat a matrix specifying which elements of the misID matrix are linked to model equations
#' @param N.par.misID a vector giving the number of parameters for each misID model (in multinomial logit space)
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits_misID<-function(DM.hab.pois,DM.hab.bern,DM.det,N.hab.pois.par,N.hab.bern.par,G.transect,Area.trans,Area.hab,Mapping,point.ind,spat.ind,grp.mean,misID,misID.mat,N.par.misID){		
  i.hurdle=1-is.null(DM.hab.bern)
  n.species=nrow(G.transect)
  n.cells=length(Area.hab)
  if(misID){
    n.misID.eq=max(misID.mat)
    MisID=vector("list",n.misID.eq)
    for(itmp in 1:n.misID.eq)MisID[[itmp]]=runif(N.par.misID[itmp],-.5,.5)
    diag.mods=diag(misID.mat)
    diag.mods=diag.mods[which(diag.mods>0)]
    if(length(diag.mods)>0){
      for(itmp in 1:length(diag.mods))MisID[[diag.mods[itmp]]][1]=MisID[[diag.mods[itmp]]][1]+2 #ensure that the highest probability is for a non-misID
    }
  }
  hab.pois=matrix(0,n.species,max(N.hab.pois.par))
  hab.bern=NULL
  tau.eta.bern=NULL
  Eta.bern=NULL
  if(i.hurdle==1){
    hab.bern=matrix(0,n.species,max(N.hab.bern.par))
    tau.eta.bern=runif(n.species,0.5,2)
    Eta.bern=matrix(rnorm(n.species*n.cells),n.species,n.cells)
  }
  Nu=matrix(0,n.species,n.cells)
  for(isp in 1:n.species){
    Nu[isp,]=log(max(G.transect[isp,])/mean(Area.trans)*exp(rnorm(length(Area.hab),0,0.1)))
  }
  Par=list(det=rnorm(ncol(DM.det),0,1),hab.pois=hab.pois,hab.bern=hab.bern,cor=ifelse(point.ind,runif(1,0,.8),0),
           Nu=Nu,Eta.pois=matrix(rnorm(n.species*n.cells),n.species,n.cells),Eta.bern=Eta.bern,
           tau.eta.pois=runif(n.species,0.5,2),tau.eta.bern=tau.eta.bern,tau.nu=runif(n.species,0.5,2),MisID=MisID)
  Par$hab.pois[,1]=log(apply(G.transect,1,'mean')/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(n.species,0,1)))
  Par$G=round(exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu))))
  for(isp in 1:n.species)Par$N[isp,]=Par$G[isp,]+rpois(n.cells,grp.mean[isp]*Par$G[isp,])
  if(spat.ind==1){
    Par$Eta.bern=0*Par$Eta.bern
    Par$Eta.pois=0*Par$Eta.pois
  }
  Par
}

#' Fill confusion array - one confusion matrix for each individual (DEPRECATED)
#' @param Confuse	An 3-dimensional array, with dimensions (# of individuals, # of rows in misID.mat, # of cols of misID.mat)
#' @param Cov	Data frame including all covariates for the misclassification model (individuals are on rows)
#' @param Beta A list where each entry is a vector giving the parameters of the misID model
#' @param n.indiv Integer giving the number of individuals 
#' @param misID.mat With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#' @param misID.formulas A formula vector providing linear model-type formulas for each positive value of misID.mat.  If the same model is used in multiple columns it is assumed that all fixed effects (except the intercept) are shared
#' @param symm	if TRUE, symmetric classification probabilities are applied (e.g. pi^12=pi^21)
#' @return A filled version of Confuse
#' @export
#' @author Paul B. Conn
get_confusion_array<-function(Confuse,Cov=NULL,Beta,n.indiv,misID.mat,misID.formulas,symm=TRUE){
  if(is.null(Cov)==1)Cov=data.frame(matrix(1,n.indiv,1))
  DM=vector("list",max(misID.mat))
  Pi=DM
  ind.mat=matrix(c(1:length(misID.mat)),nrow(misID.mat),ncol(misID.mat))
  
  for(ipar in 1:length(misID.mat)){
    if(misID.mat[ipar]==0)Pi[[ipar]]=rep(0,n.indiv)
    if(misID.mat[ipar]<0)Pi[[ipar]]=rep(1,n.indiv)
    if(misID.mat[ipar]>0){
      DM[[ipar]]=model.matrix(misID.formulas[[misID.mat[ipar]]],data=Cov)
      Pi[[ipar]]=exp(DM[[ipar]]%*%Beta[[misID.mat[ipar]]])
    }
  }
  if(symm==TRUE){
    for(iind in 1:n.indiv){
      for(icol in 1:ncol(misID.mat)){
        Confuse[iind,1,icol]=Pi[[ind.mat[1,icol]]][iind]
      }
      Confuse[iind,1,]=Confuse[iind,1,]/sum(Confuse[iind,1,])
      Pi[[ind.mat[2,3]]]=rep(1,n.indiv)
      Pi[[ind.mat[2,1]]]=(Confuse[iind,1,2]+Confuse[iind,1,2]*Pi[[ind.mat[2,2]]])/(1-Confuse[iind,1,2])
      for(icol in 1:ncol(misID.mat))Confuse[iind,2,icol]=Pi[[ind.mat[2,icol]]][iind]
      Confuse[iind,2,]=Confuse[iind,2,]/sum(Confuse[iind,2,])		
    }
  }
  else{
    for(iind in 1:n.indiv){
      for(irow in 1:nrow(misID.mat)){
        for(icol in 1:ncol(misID.mat))Confuse[iind,irow,icol]=Pi[[ind.mat[irow,icol]]][iind]
        Confuse[iind,irow,]=Confuse[iind,irow,]/sum(Confuse[iind,irow,])
      }
    }
  }
  Confuse
}

#' Fill a list with confusion matrices for each record
#' @param Cur.dat	Matrix giving data (records and covariates)  - multiple rows can be given (e.g. reflecting different observers)
#' @param stacked.names	A character vector giving column names for the data
#' @param factor.ind  An integer vector holding whehter each column of data is to be treated as numeric or factor
#' @param Levels  A list, each entry of which corresponds to a column name for factor variables and gives the possible levels of those factors
#' @param Beta A list where each entry is a vector giving the parameters of the misID model
#' @param misID.mat With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#' @param misID.models A formula vector providing linear model-type formulas for each positive value of misID.mat.  If the same model is used in multiple columns it is assumed that all fixed effects (except the intercept) are shared
#' @param misID.symm	if TRUE, symmetric classification probabilities are applied (e.g. pi^12=pi^21)
#' @return A list of confusion matrices, one for each row in Cur.dat
#' @export
#' @author Paul B. Conn
get_confusion_mat<-function(Cur.dat,Beta,misID.mat,misID.models,misID.symm=TRUE,stacked.names,factor.ind,Levels){
  Pi=vector("list",length(misID.mat))
  n.obs=nrow(Cur.dat)
  ind.mat=matrix(c(1:length(misID.mat)),nrow(misID.mat),ncol(misID.mat))
  Confuse=vector("list",n.obs)
  for(ipar in 1:length(misID.mat)){
    if(misID.mat[ipar]==0)Pi[[ipar]]=rep(0,n.obs)
    if(misID.mat[ipar]<0)Pi[[ipar]]=rep(1,n.obs)
    if(misID.mat[ipar]>0){
      DM=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=stacked.names,factor.ind=factor.ind,Det.formula=misID.models[[misID.mat[ipar]]],Levels=Levels)
      Pi[[ipar]]=exp(DM%*%Beta[[misID.mat[ipar]]])
    }
  }					
  if(misID.symm==TRUE){
    for(irow in 2:nrow(misID.mat)){
      for(icol in 1:(irow-1))Pi[[ind.mat[irow,icol]]]=rep(0,n.obs) #initialize to zero for entries set with symmetry constraint
    }
    
    for(iobs in 1:n.obs){
      Confuse[[iobs]]=matrix(0,nrow(misID.mat),ncol(misID.mat))
      #step one, calculate assignment probabilities for first row of confusion array
      for(icol in 1:ncol(misID.mat))Confuse[[iobs]][1,icol]=Pi[[ind.mat[1,icol]]][iobs]
      Confuse[[iobs]][1,]=Confuse[[iobs]][1,]/sum(Confuse[[iobs]][1,])
      #now, for remaining rows, substitute in confusion values from previous rows and calculate Pi values
      for(irow in 2:nrow(misID.mat)){
        sum.pi=0
        for(icol in 1:ncol(misID.mat))sum.pi=sum.pi+Pi[[ind.mat[irow,icol]]][iobs]
        for(icol in 1:(irow-1))Confuse[[iobs]][irow,icol]=Confuse[[iobs]][icol,irow]
        sum.Conf=sum(Confuse[[iobs]][irow,])
        for(icol in 1:(irow-1))Pi[[ind.mat[irow,icol]]][iobs]=Confuse[[iobs]][icol,irow]*sum.pi/(1-sum.Conf)
        for(icol in 1:ncol(misID.mat))Confuse[[iobs]][irow,icol]=Pi[[ind.mat[irow,icol]]][iobs]
        Confuse[[iobs]][irow,]=Confuse[[iobs]][irow,]/sum(Confuse[[iobs]][irow,])
      }
      
    }
  }
  else{
    for(iobs in 1:n.obs){
      Confuse[[iobs]]=matrix(0,dim(misID.mat))
      for(irow in 1:nrow(misID.mat)){
        for(icol in 1:ncol(misID.mat))Confuse[[iobs]][irow,icol]=Pi[[ind.mat[irow,icol]]][iobs]
        Confuse[[iobs]][irow,]=Confuse[[iobs]][irow,]/sum(Confuse[[iobs]][irow,])
      }
    }
  }
  Confuse					
}



#' compute the first derivative of log_lambda likelihood component for Langevin-Hastings
#' @param Mu	expected value for all cells
#' @param Nu	current observed valus (all cells)
#' @param Sampled Vector giving the cell identities for all sampled cells
#' @param Area Proportional area of each sampled cell that is covered by one or more transects
#' @param N	    number of groups in each transect
#' @param var.nu	variance of the overdispersion process
#' @return a gradient value
#' @export
#' @keywords gradient, Langevin-Hastings
#' @author Paul B. Conn
log_lambda_gradient<-function(Mu,Nu,Sampled,Area,N,var.nu){
  Grad=(Mu[Sampled]-Nu[Sampled])/var.nu+N-Area*exp(Nu[Sampled])
  Grad 
}

#' compute the likelihood for nu parameters
#' @param Log.lambda 	Log of poisson intensities for total areas sampled in each sampled strata
#' @param DM	the design matrix
#' @param Beta	linear predictor parameters for the log of abundance intensity	
#' @param Eta	a vector of spatial random effects
#' @param SD	standard deviation of the overdispersion process
#' @param N		a vector giving the current iteration's number of groups in the area 
#' @param Sampled Index for which cells were actually sampled
#' @param Area	Total area sampled in each sampled cell
#' @return the log likelihood associated with the data and the current set of parameters 
#' @export
#' @keywords log likelihood
#' @author Paul B. Conn
log_lambda_log_likelihood<-function(Log.lambda,DM,Beta,Eta=0,SD,N,Sampled,Area){
  Pred.log.lam=(DM%*%Beta+Eta)[Sampled]	
  logL=sum(dnorm(Log.lambda,Pred.log.lam,SD,log=1)) #normal component
  logL=logL+sum(N*Log.lambda-Area*exp(Log.lambda))
  return(logL)
} 	

#' SIMULATE AN ICAR PROCESS 
#' @param Q Precision matrix for the ICAR process
#' @return Spatial random effects
#' @export 
#' @keywords ICAR, simulation
#' @author Devin Johnson
rrw <- function(Q){
  v <- eigen(Q, TRUE)
  val.inv <- sqrt(ifelse(v$values>sqrt(.Machine$double.eps), 1/v$values, 0))
  P <- v$vectors
  sim <- P%*%diag(val.inv)%*%rnorm(dim(Q)[1], 0, 1)
  X <- rep(1,length(sim))
  if(sum(val.inv==0)==2) X <- cbind(X, 1:length(sim))
  sim <- sim-X%*%solve(crossprod(X), crossprod(X,sim))
  return(sim)
}

#' Produce an adjacency matrix for a vector
#' @param x length of vector
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn
linear_adj <- function(x){
  Adj1=matrix(0,x,x)
  Adj2=matrix(0,x,x)
  diag.min.1=diag(x-1)
  Adj1[2:x,1:(x-1)]=diag.min.1
  Adj2[1:(x-1),2:x]=diag.min.1
  Adj=Adj1+Adj2
  Adj
}	

#' Produce an adjacency matrix for a square grid
#' @param x number of cells on side of grid
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn
square_adj <- function(x){
  Ind=matrix(c(1:x^2),x,x)
  Adj=matrix(0,x^2,x^2)
  for(i in 1:x){
    for(j in 1:x){
      if(i==1 & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+x]=1
        Adj[Ind[i,j],Ind[i,j]+x+1]=1
      }
      if(i==1 & j>1 & j<x){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+x]=1
        Adj[Ind[i,j],Ind[i,j]-x]=1
        Adj[Ind[i,j],Ind[i,j]+x+1]=1
        Adj[Ind[i,j],Ind[i,j]-x+1]=1
      }
      if(i==1 & j==x){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-x]=1	
        Adj[Ind[i,j],Ind[i,j]-x+1]=1
      }
      if(i>1 & i<x & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+x]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+x-1]=1
        Adj[Ind[i,j],Ind[i,j]+x+1]=1
      }
      if(i>1 & i<x & j>1 & j<x){
        cur.nums=c(Ind[i,j]-x-1,Ind[i,j]-x,Ind[i,j]-x+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+x-1,Ind[i,j]+x,Ind[i,j]+x+1)
        Adj[Ind[i,j],cur.nums]=1
      }
      if(i>1 & i<x & j==x){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-x]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1	
        Adj[Ind[i,j],Ind[i,j]-x-1]=1
        Adj[Ind[i,j],Ind[i,j]-x+1]=1
        
      }
      if(i==x & j==1){
        Adj[Ind[i,j],Ind[i,j]+x]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1	
        Adj[Ind[i,j],Ind[i,j]+x-1]=1
      }
      if(i==x & j>1 & j<x){
        Adj[Ind[i,j],Ind[i,j]+x]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1								
        Adj[Ind[i,j],Ind[i,j]-x]=1
        Adj[Ind[i,j],Ind[i,j]+x-1]=1
        Adj[Ind[i,j],Ind[i,j]-x-1]=1
      }
      if(i==x & j==x){
        Adj[Ind[i,j],Ind[i,j]-1]=1								
        Adj[Ind[i,j],Ind[i,j]-x]=1
        Adj[Ind[i,j],Ind[i,j]-x-1]=1
      }				
    }
  }
  return(Adj)
}

#' Produce an RW1 adjacency matrix for a rectangular grid for use with areal spatial models (queens move)
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj <- function(x,y,byrow=FALSE){
  Ind=matrix(c(1:(x*y)),y,x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,x*y,x*y)
  for(i in 1:n.row){
    for(j in 1:n.col){
      if(i==1 & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i==1 & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i==1 & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i>1 & i<n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i>1 & i<n.row & j>1 & j<n.col){
        cur.nums=c(Ind[i,j]-n.row-1,Ind[i,j]-n.row,Ind[i,j]-n.row+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+n.row-1,Ind[i,j]+n.row,Ind[i,j]+n.row+1)
        Adj[Ind[i,j],cur.nums]=1
      }
      if(i>1 & i<n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
        
      }
      if(i==n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
      }
      if(i==n.row & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
      if(i==n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
    }
  }
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

#' Produce an RW2 Adjacency matrix for a rectangular grid for use with areal spatial models.
#' This formulation uses cofficients inspired by a thin plate spline, as described in Rue & Held, section 3.4.2
#' Here I'm outputting an adjacency matrix of 'neighbor weights' which makes Q construction for regular latices
#' easy to do when not trying to make inference about all cells (i.e., one can just
#' eliminate rows and columns associated with cells one isn't interested in and set Q=-Adj+Diag(sum(Adj)) 
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj_RW2 <- function(x,y,byrow=FALSE){
  cur.x=x+4  #make calculations on a larger grid and then cut off rows/columns at end
  cur.y=y+4
  Ind=matrix(c(1:(cur.x*cur.y)),cur.y,cur.x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,cur.x*cur.y,cur.x*cur.y)
  for(i in 3:(n.row-2)){
    for(j in 3:(n.col-2)){
      #kings move
      Adj[Ind[i,j],Ind[i,j]+1]=8
      Adj[Ind[i,j],Ind[i,j]+n.row]=8
      Adj[Ind[i,j],Ind[i,j]-n.row]=8
      Adj[Ind[i,j],Ind[i,j]-1]=8
      #bishops move        
      Adj[Ind[i,j],Ind[i,j]+n.row-1]=-2
      Adj[Ind[i,j],Ind[i,j]+n.row+1]=-2
      Adj[Ind[i,j],Ind[i,j]-n.row-1]=-2
      Adj[Ind[i,j],Ind[i,j]-n.row+1]=-2
      #kings move + 1
      Adj[Ind[i,j],Ind[i,j]+2]=-1
      Adj[Ind[i,j],Ind[i,j]+2*n.row]=-1  
      Adj[Ind[i,j],Ind[i,j]-2]=-1
      Adj[Ind[i,j],Ind[i,j]-2*n.row]=-1  
    }
  }
  #compile list of cells that need to be removed
  I.rem=matrix(0,n.row,n.col)
  I.rem[c(1,2,n.row-1,n.row),]=1
  I.rem[,c(1,2,n.col-1,n.col)]=1
  Adj=Adj[-which(I.rem==1),-which(I.rem==1)]
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

#' estimate optimal 'a' parameter for linex loss function
#' @param Pred.G  Predicted group abundance
#' @param Obs.G	Observed group abundance
#' @param min.a Minimum value for linex 'a' parameter
#' @param max.a Maximum value for linex 'a' parameter
#' @return The optimal tuning parameter for linex loss function as determined by minimum sum of squares 
#' @export
#' @keywords linex
#' @author Paul B. Conn
calc_linex_a<-function(Pred.G,Obs.G,min.a=0.00001,max.a=1.0){
  Y=apply(Obs.G,2,mean)
  linex_ssq<-function(a,X,Y){
    Theta=exp(-a*X)
    Theta=-1/a*log(apply(Theta,2,'mean'))
    return(sum((Y-Theta)^2))
  }
  a=optimize(f=linex_ssq,interval=c(min.a,max.a),X=Pred.G,Y=Y)
  a
} 	

#' plot 'observed' versus predicted values for abundance of each species at each transect
#' @param Out  Output list from "mcmc_ds.R" 	
#' @return NULL 
#' @export
#' @keywords diagnostics, plot
#' @author Paul B. Conn
plot_obs_pred<-function(Out){
  n.species=dim(Out$Pred.N)[1]
  par(mfrow=c(n.species,1))
  for(isp in 1:n.species){
    a.linex=calc_linex_a(Out$Pred.N[isp,,],Out$Obs.N[isp,,])$minimum
    max.x=max(c(apply(Out$Obs.N[isp,,],2,'mean'),apply(Out$Pred.N[isp,,],2,'mean')))
    plot(apply(Out$Obs.N[isp,,],2,'mean'),apply(Out$Pred.N[isp,,],2,'mean'),pch=1,xlim=c(0,max.x),ylim=c(0,max.x),xlab="Observed",ylab="Predicted")
    points(apply(Out$Obs.N[isp,,],2,'mean'),apply(Out$Pred.N[isp,,],2,'median'),pch=2)
    Theta=exp(-a.linex*Out$Pred.N[isp,,])
    Theta=-1/a.linex*log(apply(Theta,2,'mean'))
    points(apply(Out$Obs.N[isp,,],2,'mean'),Theta,pch=3)
    abline(a=0,b=1)
    legend(max.x*.1,max.x*.8,c("Mean","Median","Linex"),pch=c(1,2,3))
  }
}

#' calculate parameter estimates and confidence intervals for various loss functions
#' @param Out  Output list from "mcmc_ds.R" 	
#' @return summary.N  list vector, with the first list index indicating species
#' @export
#' @keywords summary
#' @author Paul B. Conn
summary_N<-function(Out){
  n.species=dim(Out$Pred.N)[1]
  summary.N=vector('list',n.species)
  for(isp in 1:n.species){
    a.linex=calc_linex_a(Out$Pred.N[isp,,],Out$Obs.N[isp,,])$minimum
    Theta=exp(-a.linex*Out$Post$N[isp,,])
    Theta=-1/a.linex*log(apply(Theta,2,'mean'))
    summary.N[[isp]]=list(mean=sum(apply(Out$Post$N[isp,,],2,'mean')),median=sum(apply(Out$Post$N[isp,,],2,'median')),linex=sum(Theta))
  }  
  summary.N
}

#' Mrds probit detection and related functions
#'
#' For independent observers, probit.fct computes observer-specific detection functions,
#' conditional detection functions, delta dependence function, duplicate detection function (seen by both),
#' and pooled detection function (seen by at least one).
#'
#' The vectors of covariate values can be of different lengths because expand.grid is used to create a
#' dataframe of all unique combinations of the distances and covariate values and the detection and related
#' values are computed for each combination.  The covariate vector observer=1:2 is automatically included.
#' The folowing is too long for the examples section:
#' test=probit.fct(0:10,~distance,c(1,-.15),.8,size=1:3)
#' par(mfrow=c(1,2))
#' with(test[test$observer==1,],
#' {plot(distance,p,ylim=c(0,1),xlab="Distance",ylab="Detection probability")
#' points(distance,pc,pch=2)
#' points(distance,dup,pch=3)
#' points(distance,pool,pch=4)
#' legend(1,.2,legend=c("Detection","Conditional detection","Duplicate detection","Pooled detection"),pch=1:4,bty="n")
#' plot(distance,delta,xlab="Distance",ylab="Dependence")
#' })
#' @param x vector of perpendicular distances
#' @param formula linear probit formula for detection using distance and other covariates
#' @param beta parameter values
#' @param rho maximum correlation at largest distance
#' @param ... any number of named vectors of covariates used in the formula
#' @return dat dataframe with distance, observer, any covariates specified in ... and detection probability p,
#' conditional detection probability pc, dupiicate detection dup, pooled detection pool and
#' dependence pc/p=delta.
#' @export
#' @author Jeff Laake
probit.fct=function(x,formula,beta,rho,...)
{
  #require(mvtnorm)
  #  Create dataframe and apply formula to get design matrix
  dat=expand.grid(distance=x,observer=1:2,...)
  xmat=model.matrix(formula,dat)
  #  Make sure length of beta matches number of columns of design matrix
  if(ncol(xmat)!=length(beta))stop("Mismatch between beta and formula")
  #  Compute XB and partition for 2 observers
  xbeta=xmat%*%beta
  xbeta1=xbeta[dat$observer==1]
  xbeta2=xbeta[dat$observer==2]
  #  Compute rho values
  distance=dat$distance[dat$observer==1]
  rhox=rho*distance/max(distance)
  #  Compute detection observer-specific p1,p2 and duplicate p3
  p1=pnorm(xbeta1,0,1)
  p2=pnorm(xbeta2,0,1)
  p3=apply(cbind(xbeta1,xbeta2,rhox),1,function(x)
    pmvnorm(lower=c(-x[1],-x[2]),corr=matrix(c(1,x[3],x[3],1),ncol=2,nrow=2)))
  #  Compute conditional detection prob
  p1c2=p3/p2
  p2c1=p3/p1
  #  Store values in dataframe
  dat$p[dat$observer==1]=p1
  dat$p[dat$observer==2]=p2
  dat$pc[dat$observer==1]=p1c2
  dat$pc[dat$observer==2]=p2c1
  dat$dup[dat$observer==1]=p3
  dat$dup[dat$observer==2]=p3
  dat$pool[dat$observer==1]=p1+p2-p3
  dat$pool[dat$observer==2]=p1+p2-p3
  dat$delta=dat$pc/dat$p
  return(dat)
}

#' function to convert HierarchicalDS MCMC list vector (used in estimation) into an mcmc object (cf. coda package) 
#' @param MCMC list vector providing MCMC samples for each parameter type 
#' @param N.hab.pois.par see help for mcmc_ds.R
#' @param N.hab.bern.par see help for mcmc_ds.R
#' @param Cov.par.n see help for mcmc_ds.R
#' @param Hab.pois.names see help for mcmc_ds.R
#' @param Hab.bern.names see help for mcmc_ds.R
#' @param Cov.names see help for mcmc_ds.R
#' @param Det.names see help for mcmc_ds.R
#' @param MisID.names see help for mcmc_ds.R
#' @param N.par.misID see help for mcmc_ds.R
#' @param misID.mat see help for mcmc_ds.R
#' @param misID see help for mcmc_ds.R
#' @param fix.tau.nu see help for mcmc_ds.R
#' @param spat.ind see help for mcmc_ds.R
#' @param point.ind see help for mcmc_ds.R
#' @export
#' @keywords MCMC, coda
#' @author Paul B. Conn
convert.HDS.to.mcmc<-function(MCMC,N.hab.pois.par,N.hab.bern.par,Cov.par.n,Hab.pois.names,Hab.bern.names,Det.names,Cov.names,MisID.names,N.par.misID=NULL,misID.mat=NULL,fix.tau.nu=FALSE,misID=TRUE,spat.ind=TRUE,point.ind=TRUE){
  #require(coda)
  if(misID==TRUE & (is.null(N.par.misID)|is.null(misID.mat)))cat("\n Error: must provide N.par.misID and misID.mat whenever misID=TRUE \n")
  i.ZIP=!is.na(N.hab.bern.par)[1]
  n.species=nrow(MCMC$Hab.pois)
  n.iter=length(MCMC$Hab.pois[1,,1])
  n.col=n.species*2+sum(N.hab.pois.par)+ncol(MCMC$Det)+point.ind+(1-spat.ind)*n.species+(1-fix.tau.nu)*n.species+sum(Cov.par.n)*n.species+misID*sum(N.par.misID)
  if(i.ZIP)n.col=n.col+sum(N.hab.bern.par)+(1-spat.ind)*n.species #for ZIP model
  n.cells=dim(MCMC$G)[3]
  Mat=matrix(0,n.iter,n.col)
  Mat[,1:n.species]=t(MCMC$N.tot)
  counter=n.species
  col.names=paste("Abund.sp",c(1:n.species),sep='')
  for(isp in 1:n.species){
    Mat[,counter+isp]=rowSums(as.matrix(MCMC$G[isp,,],nrow=n.iter,ncol=n.cells)) #total abundance of groups
    col.names=c(col.names,paste("Groups.sp",isp,sep=''))
  }
  counter=counter+n.species
  for(isp in 1:n.species){  #habitat parameters
    Mat[,(counter+1):(counter+N.hab.pois.par[isp])]=MCMC$Hab.pois[isp,,1:N.hab.pois.par[isp]]
    col.names=c(col.names,paste("Hab.pois.sp",isp,Hab.pois.names[[isp]],sep=''))
    counter=counter+sum(N.hab.pois.par[isp])
  }
  if(i.ZIP){
    for(isp in 1:n.species){  #habitat parameters
      Mat[,(counter+1):(counter+N.hab.bern.par[isp])]=MCMC$Hab.bern[isp,,1:N.hab.bern.par[isp]]
      col.names=c(col.names,paste("Hab.bern.sp",isp,Hab.bern.names[[isp]],sep=''))
      counter=counter+sum(N.hab.bern.par[isp])
    }   
  }
  Mat[,(counter+1):(counter+ncol(MCMC$Det))]=as.matrix(MCMC$Det)
  col.names=c(col.names,paste("Det.",Det.names,sep=''))
  counter=counter+ncol(MCMC$Det)
  if(point.ind==TRUE){
    Mat[,counter+1]=MCMC$cor
    col.names=c(col.names,"rho")
    counter=counter+1
  }
  if(spat.ind==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta.pois)
    col.names=c(col.names,paste("tau.eta.pois.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(spat.ind==FALSE & i.ZIP){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta.bern)
    col.names=c(col.names,paste("tau.eta.bern.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(fix.tau.nu==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.nu)
    col.names=c(col.names,paste("tau.nu.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(is.null(Cov.par.n)==FALSE){
    max.par=max(Cov.par.n)
    for(isp in 1:n.species){
      for(ipar in 1:length(Cov.par.n)){
        Mat[,(counter+1):(counter+Cov.par.n[ipar])]=MCMC$Cov.par[isp,,((ipar-1)*max.par+1):((ipar-1)*max.par+Cov.par.n[ipar])]
        counter=counter+Cov.par.n[ipar]
        col.names=c(col.names,paste("Cov.sp",isp,".",Cov.names[[ipar]],sep=''))
      }
    }
  }
  if(misID==TRUE){
    for(imod in 1:max(misID.mat)){
      Mat[,(counter+1):(counter+N.par.misID[imod])]=MCMC$MisID[[imod]]
      counter=counter+N.par.misID[imod]
      col.names=c(col.names,paste("misID.mod",imod,".",MisID.names[[imod]],sep=''))
    }
  }	
  colnames(Mat)=col.names
  Mat=mcmc(Mat)
  Mat
}


#' function to export posterior summaries from an mcmc object to a table
#' @aliases table.mcmc
#' @S3method table mcmc
#' @method table mcmc
#' @param MCMC An mcmc object with columns referencing different parameter types (column names are used for plotting labels)
#' @param file A file name to ouput to (including path); if null (default), outputs to screen
#' @param type What type of table to produce (either "csv" or "tex")
#' @param a Value to use for credible intervals.  For example, alpha=0.05 results in 95\% credible intervals
#' @export
#' @keywords MCMC, table
#' @author Paul B. Conn
table.mcmc<-function(MCMC,file=NULL,type="csv",a=0.05){
  #require(xtable)
  Out.tab=data.frame(matrix(0,ncol(MCMC),5))
  colnames(Out.tab)=c("Parameter","Mean","Median","Lower","Upper")
  MCMC=as.matrix(MCMC)
  Out.tab[,1]=colnames(MCMC)
  Out.tab[,2]=colMeans(MCMC)
  Out.tab[,3]=apply(MCMC,2,'median')
  Out.tab[,4]=apply(MCMC,2,'quantile',a/2)
  Out.tab[,5]=apply(MCMC,2,'quantile',1-a/2)
  if(is.null(file))print(Out.tab)
  else{
    if(type=="csv")write.csv(Out.tab,file=file)
    if(type=="tex"){
      Out.tab=xtable(Out.tab)
      print(Out.tab,file=file)
    }
    if(type!="csv" & type!="tex")cat("\n Error: unknown table type.  No table was printed to file.")
  }
}

#' function to calculate posterior predictive loss given the output object from hierarchicalDS
#' @param Out Output object from running hierarchicalDS
#' @param burnin Any additional #'s of values from beginning of chain to discard before calculating PPL statistic (default is 0)
#' @return A matrix with posterior variance (P), sums of squares (G) for the posterior mean and median predictions (compared to Observations), and total posterior loss (D)
#' @export
#' @keywords Posterior predictive loss
#' @author Paul B. Conn
post_loss<-function(Out,burnin=0){
  dims.Pred=dim(Out$Pred.det)
  median.Pred=array(0,dim=dims.Pred[2:4])
  mean.Pred=median.Pred
  var.Pred=mean.Pred
  for(itrans in 1:dims.Pred[2]){
    for(isp1 in 1:dims.Pred[3]){
      for(isp2 in 1:dims.Pred[4]){
        median.Pred[itrans,isp1,isp2]=median(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,isp1,isp2])
        mean.Pred[itrans,isp1,isp2]=mean(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,isp1,isp2])
        var.Pred[itrans,isp1,isp2]=var(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,isp1,isp2])
      }
    }  
  }
  sum.sq.mean=sum((Out$Obs.det-mean.Pred)^2)
  sum.sq.median=sum((Out$Obs.det-median.Pred)^2)
  Loss=matrix(0,2,3)
  colnames(Loss)=c("P","G","D")
  rownames(Loss)=c("mean","median")
  Loss[,1]=sum(var.Pred)
  Loss[1,2]=sum.sq.mean
  Loss[2,2]=sum.sq.median
  Loss[,3]=rowSums(Loss[1:2,1:2])
  Loss
}

#' function to plot a map of abundance.  this was developed for spatio-temporal models in mind
#' @param cur.t time step to plot
#' @param N A vector of values to plot (e.g. abundance) 
#' @param Grid A list of SpatialPolygonsDataFrame (one for each time step) - holding survey unit spatial information 
#' @param highlight If provided, the rows of Grid[[cur.t]] to specially highlight
#' @return A ggplot2 object
#' @export
#' @keywords Abundance map
#' @author Paul B. Conn
plot_N_map<-function(cur.t,N,Grid,highlight=NULL){
  #require(rgeos)
  Tmp=Grid[[1]]
  if(is.null(highlight)==FALSE){
    midpoints=data.frame(gCentroid(Tmp[highlight,],byid=TRUE))
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N[,cur.t]
  Cur.df=cbind(data.frame(gCentroid(Tmp,byid=TRUE)),Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  p1=ggplot(Cur.df)+aes(x=Easting,y=Northing,fill=Abundance)+geom_raster()+tmp.theme
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
  }
  p1
}

#' MCMC output from running example in Hierarchical DS 
#' 
#' @name sim_out 
#' @docType data 
#' @author Paul Conn \email{paul.conn@@noaa.gov} 
#' @keywords data 
NULL