#' function to simulate distance sampling data from simple model with increasing abundance
#' intensity, no assumed spatial structure, and point independence.  If no parameters are given, uses internally defined values.
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @param Observers A (2 x S) matrix giving the observer identity for each transect
#' @param misID If TRUE (default), misidentification is assumed
#' @param X.site model.matrix for habitat covariates (defaults to a linear, quadratic effects of transect # on log scale)
#' @param n.species Number of species to simulate (current max is 2) (default is 2)
#' @param Beta.hab A (# of species X # parameters) matrix giving parameters for habitat-abundance relationship (default is linear increase for species 1, quadratic for species 2)
#' @param detect.model	A formula for the detection model.  Default is ~Observer+Distance+Group+Species (formula should consist of these key words)
#' @param Beta.det A vector giving parameters for the detection function; # of parameters must match model.matrix()!
#' @param dist.cont If TRUE, uses continuous distances on (0,1).  If FALSE, uses discrete distance classes
#' @param n.bins If dist.cont=FALSE, how many bins to use for distances.  Default is 5.
#' @param cor.par Correlation at maximum distance.  Default is 0.5
#' @param Grp.par	A vector with an entry for each species giving parameters for group size (assumed zero-truncated Poisson). Default is 3 and 1, corresponding to mean group sizes of 4 and 2 for each species
#' @param misID.mat With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#' @param misID.par A list, each element of which gives the parameters associated with each entry in misID.models 
#' @param misID.models A formula vector providing linear model-type formulas for each positive value of misID.mat.  If the same model is used in multiple columns it is assumed that all fixed effects (except the intercept) are shared
#' @param misID.symm If TRUE, the constraint pi^{i|j}=pi^{j|i} is implemented; in this case, entries for pi^{j|i} are all assumed to be = pi^{i|j} (default is TRUE)
#' @return a distance sampling dataset
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_data_overd<-function(S,Observers,misID=TRUE,X.site=NULL,n.species=2,Beta.hab=NULL,Beta.det=NULL,detect.model=~Observer+Distance+Group,dist.cont=FALSE,n.bins=5,cor.par=0.5,Grp.par=c(3,1),misID.models=NULL,misID.par=NULL,misID.mat=NULL,misID.symm=TRUE){
	require(mvtnorm)
	
	if(n.species>2)cat("\n Error: current max species is 2 \n")
	if(n.species==1 & misID==TRUE){
		cat("\n n.speces=1 so misID set to FALSE \n")
		misID=FALSE
	}
	if((n.bins!=5 | dist.cont==TRUE) & (is.null(Beta.det)==TRUE))cat("\n Error: if using continuous distances or a non-default distance bin #, you must input Beta.det \n")
	#process parameters
	if(is.null(X.site)==TRUE)X.site=cbind(rep(1,S),log(c(1:S)/S),(log(c(1:S)/S))^2) #covariate on abundance intensity, sp 1
	if(is.null(Beta.hab)==TRUE){
		Beta.hab=matrix(0,n.species,3)
		#Beta.hab[1,]=c(log(40),1,0) 
		Beta.hab[1,]=c(log(100),1,0) 
		if(n.species==2)Beta.hab[2,]=c(log(10),-2,-1)
		Beta.hab[1,]=c(log(40),0,0) 
		if(n.species==2)Beta.hab[2,]=c(log(10),0,0)
	}
	
	#detection parameters
	if(dist.cont==FALSE)Levels=list(Observer=sort(unique(c(Observers))),Distance=as.factor(c(1:n.bins)),Species=as.factor(1:n.species))
	else Levels=list(Observer=unique(c(Observers)),Species=as.factor(1:n.species))
	factor.ind=list(Observer=TRUE,Distance=(dist.cont==FALSE),Group=FALSE,Species=TRUE)
	
	#if(is.null(Beta.det)==TRUE)Beta.det=c(10,-.2,-.4,-.6,-.9,-1.1,-1.3,.1,.3)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
	if(is.null(Beta.det)==TRUE)Beta.det=c(1.2,-.2,-.4,-.6,-.9,-1.1,-1.3,.1)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size
	
  N1=c(0,0,0,0,0,0,1,0,0,0,2,0,0,14,4,3)*3
  N2=c(0,5,0,7,0,0,0,25,3,1,12,3,0,0,0,0)*3
  #N1=round(exp(X.site%*%Beta.hab[1,]))
	#N2=N1*0
	#if(n.species==2)N2=round(exp(X.site%*%Beta.hab[2,]))
	cat(paste("\n True G, sp 1 = ",N1,'\n\n'))
	cat(paste("\n True G.tot, sp 1= ",sum(N1),'\n'))
	if(n.species==2){
		cat(paste("\n True G, sp 2 = ",N2,'\n\n'))
		cat(paste("\n True G.tot, sp 2= ",sum(N2),'\n'))
	}
	Bin.length=c(1,1.2,1.5,2,3)
	
	Dat=matrix(0,sum(N1)+sum(N2),8)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size
	
	#initialize confusion matrix, etc.
	if(misID==1){
		if(is.null(misID.mat)){
			misID.mat=matrix(0,2,3)
			misID.mat[1,]=c(1,-1,2)
			#Confusion[2,]=c(-1,1,2)
			misID.mat[2,]=c(-1,3,-1)
		}
		if(is.null(misID.models)==TRUE)misID.models=c(~1,~1,~1)
		#Conf.model=c(~Observer+Group+Distance+Species,~Observer)
		if(is.null(misID.par)==TRUE){misID.par=list(2,1,3)
			#Conf.par=list(c(2,.2,.4,.3,-.1,-.2,-.4,-.8,.5),c(-1,-.2,.2)) #parameters for getting an 'unknown'
		}
	}
	
	pl=1
	for(i in 1:S){
		cur.Observers=Observers[,i]
		n.observers=2-is.na(cur.Observers[2])
		if(N1[i]>0){
			for(j in 1:N1[i]){
				if(dist.cont==TRUE)dist=runif(1)
				else dist=sample(c(1:n.bins),1,prob=Bin.length)
				grp.size=rpois(1,Grp.par[1])+1
				Dat1=matrix(c(cur.Observers[1],dist,grp.size,1),1,4)
				if(n.observers==2){
          Dat2=Dat1
				  Dat2[1]=cur.Observers[2]
          X2=get_mod_matrix(Cur.dat=Dat2,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				}
				X1=get_mod_matrix(Cur.dat=Dat1,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=dist
				Dat[pl,7]=grp.size
				Dat[pl,8]=1
				if(dist.cont==FALSE)cur.cor=(dist-1)/(n.bins-1)*cor.par
				else cur.cor=dist*cor.par
				mu1=X1%*%Beta.det
        if(n.observers==2){
				  mu2=X2%*%Beta.det
				  Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
        }
        else Dat[pl,4:5]=c(rnorm(1,mu1,1),NA)
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		if(N2[i]>0){
			for(j in 1:N2[i]){
				if(dist.cont==TRUE)dist=runif(1)
				else dist=sample(c(1:n.bins),1,prob=Bin.length)
				grp.size=rpois(1,Grp.par[2])+1
				Dat1=matrix(c(cur.Observers[1],dist,grp.size,2),1,4)
				X1=get_mod_matrix(Cur.dat=Dat1,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				if(n.observers==2){
				  Dat2=Dat1
				  Dat2[1]=cur.Observers[2]
				  X2=get_mod_matrix(Cur.dat=Dat2,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				}
        Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=dist
				Dat[pl,7]=grp.size
				Dat[pl,8]=2
				mu1=X1%*%Beta.det
				if(dist.cont==FALSE)cur.cor=(dist-1)/(n.bins-1)*cor.par
				else cur.cor=dist*cor.par
        if(n.observers==2){
          mu2=X2%*%Beta.det
          Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))          
        }
        else Dat[pl,4:5]=c(rnorm(1,mu1,1),NA)         
        Dat[pl,4:5]=c(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		
	}
	cat(paste("Total N, sp 1 = ",sum(Dat[Dat[,8]==1,7])))
	if(n.species==2)cat(paste("Total N, sp 2 = ",sum(Dat[Dat[,8]==2,7])))
	Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),] #get rid of animals never observed
	
	#put things in "Jay's" format
	Dat2=matrix(0,2*nrow(Dat),7)
	ipl=1
	for(irecord in 1:nrow(Dat)){
		Dat2[ipl,1]=Dat[irecord,1]
		Dat2[ipl+1,1]=Dat[irecord,1]
		Dat2[ipl,3]=Dat[irecord,2]
		Dat2[ipl+1,3]=Dat[irecord,3]
		Dat2[ipl,4]=Dat[irecord,4]
		Dat2[ipl+1,4]=Dat[irecord,5]
		Dat2[ipl,5]=Dat[irecord,6]
		Dat2[ipl+1,5]=Dat[irecord,6]
		Dat2[ipl,6]=Dat[irecord,7]
		Dat2[ipl+1,6]=Dat[irecord,7]
		Dat2[ipl,7]=Dat[irecord,8]
		Dat2[ipl+1,7]=Dat[irecord,8]
		Dat2[ipl,2]=irecord  #match number
		Dat2[ipl+1,2]=irecord
		ipl=ipl+2
	}
	Dat2=as.data.frame(Dat2)
	colnames(Dat2)=c("Transect","Match","Observer","Obs","Distance","Group","Species")
	Dat2[,"Observer"]=as.factor(Dat2[,"Observer"])
	Dat2[,"Distance"]=as.factor(Dat2[,"Distance"])
  
  #get rid of NA records for transects where there was only one observer
  if(sum(is.na(Dat2[,4]))>0){
    which.NA=which(is.na(Dat2[,4]))
    Dat2=Dat2[-which.NA,]
  }
	
	Dat=Dat2
	True.species=Dat[,"Species"]
	n.indiv=nrow(Dat)
	
	subset_conf_array2<-function(Confuse,Species,n.indiv){
		Cur=matrix(0,n.indiv,dim(Confuse)[3])
		for(iind in 1:n.indiv)Cur[iind,]=Confuse[iind,Species[iind],]
		Cur
	}
	if(misID==TRUE){
	  Confuse=array(0,dim=c(n.indiv,dim(misID.mat)))
	  Confuse=get_confusion_array(Confuse,Cov=NULL,Beta=as.matrix(misID.par),n.indiv=n.indiv,misID.mat=misID.mat,misID.formulas=misID.models,symm=misID.symm)
		# Now, put in partial observation process
		Probs=subset_conf_array2(Confuse=Confuse,Species=Dat[,"Species"],n.indiv=n.indiv)		
		get_samp<-function(prob)sample(c(1:length(prob)),1,prob=prob)
		Dat[,"Species"]=apply(Probs,1,'get_samp')	
	}
	Dat[,"Species"]=as.integer(Dat[,"Species"])
	
	Dat=cbind(Dat[,1:4],Dat[,"Species"],Dat[,5:6])
	colnames(Dat)[5]="Species"
	G.true=N1
    if(n.species==2)G.true=cbind(G.true,N2)
	
	Out=list(Dat=Dat,G.true=G.true,True.species=True.species)
}

