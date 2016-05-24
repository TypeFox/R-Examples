#' function to simulate double observer spatial distance sampling data subject to possible zero inflation and species misidentification
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @param Observers A (2 x S) matrix giving the observer identity for each transect
#' @param ZIP If TRUE, simulate abundance using a zero-inflated Poisson model
#' @param misID If TRUE, assume species misidentification
#' @param tau.pois Precision of the ICAR process for Poisson abundance
#' @param tau.bern Precision of the ICAR process for zero inflation
#' @return a distance sampling dataset
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_data<-function(S,Observers,ZIP=TRUE,misID=TRUE,tau.pois=15,tau.bern=20){
  #note: currently hardwired for 2 species
	#require(mvtnorm)
	#set.seed(122412)
	
  if((sqrt(S)%%1)!=0)cat("Error: S must be a square #")
	Adj=square_adj(sqrt(S))
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q.pois=Matrix(tau.pois*Q)
  Q.bern=Matrix(tau.bern*Q)
	#simulate icar process
	Eta.pois.sp1=rrw(Q.pois)
  Eta.pois.sp2=rrw(Q.pois)
  Eta.bern.sp1=rrw(Q.bern)
  Eta.bern.sp2=rrw(Q.bern)
	#Eta=rep(0,S)
	
	SP.pois.sp1=matrix(Eta.pois.sp1,sqrt(S),sqrt(S))
	SP.pois.sp2=matrix(Eta.pois.sp2,sqrt(S),sqrt(S))
	SP.bern.sp1=matrix(Eta.bern.sp1,sqrt(S),sqrt(S))
	SP.bern.sp2=matrix(Eta.bern.sp2,sqrt(S),sqrt(S))
	
	#process parameters
	lambda.grp1=3
	lambda.grp2=1
	X.site1=cbind(rep(1,S),rep(log(c(1:sqrt(S)/sqrt(S))),each=sqrt(S))) #covariate on abundance intensity, sp 1
	X.site2=X.site1 #covariate on abundance intensity, sp 1
	Beta.site1=c(log(40),1) 
	Beta.site2=c(log(10),0)
  X.bern1=matrix(1,S,1)
  X.bern2=X.bern1
  Beta.bern1=-10
  Beta.bern2=-10
  if(ZIP){
    Beta.bern1=0
    Beta.bern2=0
  }
	
	#detection parameters
	n.bins=5 #n.bins=5 hardwired elsewhere
	Beta.det=c(1.6,1.4,1.2,-.8,-.6,-.4,-.2,.2,0,0,0,0)
	#Beta.det=c(10,10,10,-.6,-.3,-.2,-.2,.1,0,.3)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
												#in this version, distance pars are additive (i.e., bin 3 gets bin 2 and bin 3 effect).
	cor.par=0.5 #correlation in max age bin (linear from zero)
		
	#N=rpois(S,0.5*exp(X.site%*%Beta.site))
	N1=rbern(S,pnorm(X.bern1%*%Beta.bern1+as.vector(SP.bern.sp1)))*(rpois(S,exp(X.site1%*%Beta.site1+as.vector(SP.pois.sp1))))
	N2=rbern(S,pnorm(X.bern2%*%Beta.bern2+as.vector(SP.bern.sp2)))*(rpois(S,exp(X.site2%*%Beta.site2+as.vector(SP.pois.sp2))))
	cat(paste("\n True G, sp 1 = ",N1,'\n\n'))
	cat(paste("\n True G.tot, sp 1= ",sum(N1),'\n'))
	cat(paste("\n True G, sp 2 = ",N2,'\n\n'))
	cat(paste("\n True G.tot, sp 2= ",sum(N2),'\n'))
	
	Dat=matrix(0,sum(N1)+sum(N2),8)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size
	misID.mat=matrix(0,2,3)  # misID matrix 
	misID.mat[1,]=c(1,-1,2)  # positive numbers specify a model, 0 denotes impossible, -1 denotes obtain by subtraction
	misID.mat[2,]=c(-1,3,-1)
	misID.models=c(~1,~1,~1)
	MisID=vector("list",max(misID.mat))
	MisID[[1]]=2 #parameters for getting it right
	MisID[[2]]=1
	MisID[[3]]=3 #parameters for getting it right
	misID.symm=TRUE
	
	X=rep(0,length(Beta.det))
	pl=1
	for(i in 1:S){
		cur.Observers=Observers[,i]
		if(N1[i]>0){
			for(j in 1:N1[i]){
				X1=X
				X2=X
				X1[cur.Observers[1]]=1
				X2[cur.Observers[2]]=1
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=sample(c(1:n.bins),1)
				Dat[pl,7]=rpois(1,lambda.grp1)+1
				cur.sp=1
				Dat[pl,8]=cur.sp
				if(Dat[pl,6]>1){
					X1[4:(2+Dat[pl,6])]=1
					X2[4:(2+Dat[pl,6])]=1
				}
				X1[8]=Dat[pl,7]
				X2[8]=Dat[pl,7]
				temp=c(0,0)
				temp[cur.sp]=1
				X1[9:10]=temp
				X2[9:10]=temp
				mu1=X1%*%Beta.det
				mu2=X2%*%Beta.det
				cur.cor=(Dat[pl,6]-1)/(n.bins-1)*cor.par
				Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		if(N2[i]>0){
			for(j in 1:N2[i]){
				X1=X
				X2=X
				X1[cur.Observers[1]]=1
				X2[cur.Observers[2]]=1
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=sample(c(1:n.bins),1)
				Dat[pl,7]=rpois(1,lambda.grp2)+1
				cur.sp=2
				Dat[pl,8]=cur.sp
				if(Dat[pl,6]>1){
					X1[4:(2+Dat[pl,6])]=1
					X2[4:(2+Dat[pl,6])]=1
				}
				X1[8]=Dat[pl,7]
				X2[8]=Dat[pl,7]
				temp=c(0,0)
				temp[cur.sp]=1
				X1[9:10]=temp
				X2[9:10]=temp
				mu1=X1%*%Beta.det
				mu2=X2%*%Beta.det
				cur.cor=(Dat[pl,6]-1)/(n.bins-1)*cor.par
				Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		
	}
	cat(paste("Total N, sp 1 = ",sum(Dat[Dat[,8]==1,7])))
	cat(paste("Total N, sp 2 = ",sum(Dat[Dat[,8]==2,7])))
	#Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),]
	
	#put things in "Jay's" format
	Dat2=rbind(Dat,Dat)
	ipl=1
	for(irecord in 1:nrow(Dat)){
		Dat2[ipl,1]=Dat[irecord,1]
		Dat2[ipl+1,1]=Dat[irecord,1]
		Dat2[ipl,3]=Dat[irecord,2]
		Dat2[ipl+1,3]=Dat[irecord,3]
		Dat2[ipl,4]=Dat[irecord,4]
		Dat2[ipl+1,4]=Dat[irecord,5]
		Dat2[ipl,5]=1  #observer covariate that has no effect
		Dat2[ipl+1,5]=0
		Dat2[ipl,6]=Dat[irecord,6]
		Dat2[ipl+1,6]=Dat[irecord,6]
		Dat2[ipl,7]=Dat[irecord,7]
		Dat2[ipl+1,7]=Dat[irecord,7]
		Dat2[ipl,8]=Dat[irecord,8]
		Dat2[ipl+1,8]=Dat[irecord,8]
		Dat2[ipl,2]=irecord  #match number
		Dat2[ipl+1,2]=irecord
		ipl=ipl+2
	}
	Dat2=as.data.frame(Dat2)
	colnames(Dat2)=c("Transect","Match","Observer","Obs","Seat","Distance","Group","Species")
	Dat2[,"Observer"]=as.factor(Dat2[,"Observer"])
	Dat2[,"Distance"]=as.factor(Dat2[,"Distance"])
	Dat2[,"Seat"]=as.factor(Dat2[,"Seat"])
  #Dat2[,"Species"]=as.factor(Dat2[,"Species"])  not set to factor at this point since species is sampled below
	
	Dat=Dat2
	True.species=Dat[,"Species"]
	
  stacked.names=colnames(Dat)
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
	
	Conf=get_confusion_mat(Cur.dat=Dat,Beta=MisID,misID.mat=misID.mat,misID.models=misID.models,misID.symm=misID.symm,stacked.names=stacked.names,factor.ind=factor.ind,Levels=Levels)  				
	
	if(misID==TRUE){
		# Now, put in partial observation process
		Ind.sp1=which(Dat[,"Species"]==1)
		Ind.sp2=which(Dat[,"Species"]==2)
		Dat1=Dat[Ind.sp1,]
    Psi=matrix(0,nrow(Dat1),3)
    for(i in 1:nrow(Dat1))Psi[i,]=Conf[[Ind.sp1[i]]][1,]
    get_samp<-function(prob)sample(c(1:length(prob)),1,prob=prob)
		Cur.sp=apply(Psi,1,'get_samp')	
		Dat[Ind.sp1,"Species"]=Cur.sp
		
		Dat1=Dat[Ind.sp2,]
		Psi=matrix(0,nrow(Dat1),3)
		for(i in 1:nrow(Dat1))Psi[i,]=Conf[[Ind.sp1[i]]][2,]
		Cur.sp=apply(Psi,1,'get_samp')	
		Dat[Ind.sp2,"Species"]=Cur.sp
	}
	Dat[,"Species"]=as.integer(Dat[,"Species"])
	
	Dat=cbind(Dat[,1:4],Dat[,"Species"],Dat[,5:7])
	colnames(Dat)[5]="Species"
    G.true=cbind(N1,N2)
	
	Out=list(Dat=Dat,G.true=G.true,True.species=True.species)
}

