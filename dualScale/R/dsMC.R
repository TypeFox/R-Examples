dsMC <-
function(X,dim=NA){
  gc()
  memoryclean<-5999 #
	NR<-nrow(X) #number of respondents
	N<-ncol(X) #number of multipleChoice items
	X<-matrix(as.matrix(X),nrow=NR,ncol=N)
	K<-1 # plain DS

  itemOp<-apply(X,2,max)
  M<-sum(itemOp) # number of options
#
  if(sum(unlist(apply(X,2,tabulate))==0)>0){stop('dsHELP: You have to use dsCHECK before dsMC')}
  #if(sum(unlist(apply(X,2,tabulate))==1)>0){stop('dsHELP: There is a variable with only 1 answer')} 
  if(sum(is.na(X))>0){stop('dsHELP: You have NA values')}
  #
  dimnames(X)[[2]]<-paste("q.",1:N,sep="")
  dimnames(X)[[1]]<-paste("s.",1:NR, sep="")
  ItOpNa<-c()
  for(j in 1:N){
    ION<-paste(dimnames(X)[[2]][j], 1:itemOp[j], sep=":")
    ItOpNa<-c(ItOpNa,ION)	
    }
#	

	#
  if(is.na(dim))  {solMC<-(M-N)} else{solMC<-dim} 
  #
  F<-Matrix(0,NR,M,sparse=TRUE)
  auxiliar<-{}
  auxiliar[1]<-0
  for(j in 2:N){
    auxiliar[j]<-auxiliar[j-1]+itemOp[j-1]
  }
  F=sparseMatrix(x=1,i=c(1:NR)%*%t(matrix(1,N,1)),j=t(auxiliar+t(X)))
  auxiliar[N+1]<-auxiliar[N]+itemOp[N]
  C=t(F)%*%F
  C=C/(N-1+K)
	if(solMC>(M-N))
		stop("dsHELP: The number of dimensions have to be less than or equal to ",(M-N))
		#
  f_t<-NR*(N-1+K)
  f_c<-apply(F,2,sum) #the column sum
  f_r<-apply(F,1,sum)# the row sum
  CorrTerm<-(f_c^.5%*%t(f_c^.5))/f_t
		# 
	ola<-svd(((diag(f_c^-.5)%*%C%*%diag(f_c^-.5))-CorrTerm))
  rm(C);rm(CorrTerm);
  if (NR>memoryclean){for (i in 1:10){gc()}} ## Cleaning the memory garbarge
		#
#Table of Statistics. Step [A]
		#
	eigenval<-ola$d
	#	
control<-sum(eigenval>1/10^15)
if(control<solMC){
  maxi<-sum(eigenval>0.000001)
  stop('dsHELP: due to the eigenvalues, your data allow you to use up to ',maxi ,' components')}
singuval<-eigenval^.5
  delta<-eigenval/sum(eigenval)*100
  result<-data.frame(Component=c(1:M),Eigenvalue=round(eigenval,4), SingValue=singuval, Alpha=(1-(1-eigenval)/((N-1)*eigenval)), Delta=delta, CumDelta=cumsum(delta)/sum(delta)*100)
  #
#Options weights
		#
  w<-(f_t/diag(t(ola$u)%*%ola$u))^.5*ola$u[,1:(M-N)]
	x=diag(f_c^-.5)%*%w
	xAdj=x%*%diag((eigenval^.5)[1:(M-N)])
  rm(f_c)
  yAdj=(1/f_r)*F%*%x

  y=yAdj%*%diag(1/singuval[1:(M-N)])
#
#The correlation matrix. STEP[C]
#
  if (NR>memoryclean){for (i in 1:10){gc()}}

  cm1<-ff(initdata=0,dim=c(NR,N,(M-N)),dimorder=c(1,2,3))
  for(j in 1:(M-N)){
    q<-t(F)*xAdj[,j]
    q<-as.matrix(q)
    q<-matrix(q[(q!=0)],ncol=N,nrow=NR,byrow=TRUE)
    cm1[,,j]<-q
  }
  rm(F);rm(f_r);rm(q);
  if (NR>memoryclean){for (i in 1:10){gc()}} ## Cleaning the memory garbarge

  cm<-ff(initdata=0,dim=c(NR+1,N+1,(solMC)),dimorder=c(1,2,3))
  cm[1:NR,1:N,]<-cm1[,,1:(solMC)]
  cm[1:NR,N+1,]<-ffapply(MARGIN=c(1,3),AFUN="sum",X=cm1[1:NR,1:N,1:(solMC)], CFUN="cbind",RETURN=TRUE)##[,(1:solMC)]
  rm(cm1);
  if (NR>memoryclean){for (i in 1:10){gc()}} ## Cleaning the memory garbarge
  cm[NR+1,,]<-ffapply(MARGIN=c(2,3),AFUN="sum",X=cm[1:NR,,], CFUN="cbind",RETURN=TRUE)[,]
  if (NR>memoryclean){for (i in 1:10){gc()}} ## Cleaning the memory garbarge
		#
		#
		# STEP [D]
		#
    corM<-ff(initdata=0,dim=c(N+1,N+1,(solMC)),dimorder=c(1,2,3))
		for(i in 1:(solMC)){
				corM[,,i]<-cor(cm[1:NR,,i]) 
				}
  if (NR>memoryclean){for (i in 1:10){gc()}} ## Cleaning the memory garbarge}
  tabE1<-t(corM[,N+1,])^2
	tabE1[,N+1]<-apply(tabE1[,1:N],1,mean)
  
	tabE1ok<-tabE1
	tabE1ok[is.na(tabE1ok)]<-0 
	sumTabE1<-apply(tabE1ok,2,sum)
	labels<-paste("Item",1:(N+1))
		#
	tabE1<-round(tabE1,3)
		#
		#
	dimnames(tabE1)[[2]]<-c(paste("q.",1:N,sep=""),"Avge")
		#############
  tabE2<-array(0,dim=c(N+1,3,(solMC)))#Esta es SSj
  tabE2[,1,]<-t(t(apply(cm[1:NR,,]^2,c(2,3),sum))/eigenval[1:solMC])
  tabE2[,2,]<-t(tabE1)
  tabE2[,3,]<-tabE2[,2,]^.5
  rm(cm)

  Table<-array(0,dim=c(N,4,solMC))
  Table[,1,]<-c(1:N)
  Table[,2:4,]<-round(tabE2[1:N,,],4)
  rm(tabE2)
  #
  cramer<-array(0,dim=c(N,N))#Esta es SSj
  for(i in 1:N){
    for(j in 1:N){
      cramer[i,j]<-round(assocstats(table(X[,i],X[,j]))$cramer,2)
    }
  }

#
tipo<-"O" #for O-rdinary Dual Scaling of Multiple Choice Data
#  
		# ############OUTPUT
		
	dsMC.output<-list(
	   IniDat=X,
	   Call=match.call(),
	   CramerV=cramer,
	   eigenval=eigenval,
     SingValue=as.matrix(eigenval^.5),
	   ItONa=ItOpNa,
	   N.Comp=(solMC),
	   N.Item=N,
     N.Ss=NR,	    
     NoOpI=itemOp,
	   SubNa=dimnames(X)[[1]], 
	   tipo=tipo,
	   Tot.Op=M,
     #######
	   Proj.Op_O=round(xAdj,3),
	   Proj.Su_O=round(yAdj,3),
	   Inf_O=as.data.frame(tabE1),
	   ItemStat_O=round(Table,3),
	   Out_O=round(result[1:(M-N),],3), 
	   Rij_O=round(corM[1:N,1:N,],3), 
	   Norm.Op_O=round(x,3),
	   Norm.Su_O=round(y,3) 
)
class(dsMC.output)<-"ds"	
return(dsMC.output)}
