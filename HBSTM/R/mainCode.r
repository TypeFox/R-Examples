


##################################################################################################################
#######################################      Main functions      #################################################
##################################################################################################################


.HBSTM.creation=function(Zt,K,newGrid,reglag,seas,spatlags,hyperpriors,initialvalues,nIter,nBurn,fit,plots,posterior,save,control){
	
	#Aqui poner control de atributos de la función
	if(missing(posterior)) posterior="mean"
	controlfun(Zt=Zt,K=K,newGrid=newGrid,reglag=reglag,spatlags=spatlags,posterior)
	
	if (missing(control)) control=list()
	control  <- do.call(hbstmControl, control)

	#Control de la longitud de los parámetros espaciales

	cols=length(levels(as.factor(newGrid[,1])))
	rows=length(levels(as.factor(newGrid[,2])))
	lengthSpatPar=c(rows-1,cols-1,rep(min(rows-1,cols-1),2))
	lengthSpatPar=controlspatlags(lengthSpatPar,spatlags)
		
	#Creación de la estructura de Hyperpriors
	if(missing(hyperpriors)){
		Hyperpriors=HyperStandard.creation(reglag=reglag,Nseas=length(seas),spatlags=lengthSpatPar,S=dim(newGrid)[1])
	}else{
		Hyperpriors=Hyperpriors.control(hyperpriors,reglag=reglag,Nseas=length(seas),spatlags=lengthSpatPar)
	}	

	#Creación de la estructura de parámetros
	Parameters=Parameter.creation(reglag=reglag,seas=seas,spatlags=spatlags,S=dim(newGrid)[1],T=dim(Zt)[2])

	
	#Asignación de los valores iniciales
	if(missing(initialvalues)){
		Parameters=initialStand.creation(Parameters=Parameters,Zt=as.matrix(Zt),P=as.matrix(cbind(rep(1,dim(newGrid)[1]),as.matrix(newGrid))))
	}else{
		Parameters=initial.control(Parameters=Parameters,initialvalues=initialvalues)
	}

	if(posterior=="mean"){
		n3=2
		name=c("mean","sd")
	}else{
		n3=3
		name=c("median","LCI.95","UCI.95")
	}

	HBSTM=new(Class="HBSTM",Parameters=Parameters,Hyperpriors=Hyperpriors,newGrid=as.matrix(newGrid),seed=control$seed,iterations=0,K=K,Zt=as.matrix(Zt),fitted=array(NA,dim=c(dim(newGrid)[1],dim(Zt)[2],n3),dimnames=list(paste("Spatial",1:dim(newGrid)[1]),paste("Temp",1:dim(Zt)[2]),name)),residuals=matrix(NA,nrow=dim(Zt)[1],ncol=dim(Zt)[2]),MCMCsamp=list(NULL),MCMCclass="empty")

	if(fit){
		if(missing(nIter)) nIter=1000
		if(missing(nBurn)) nBurn=as.integer(nIter*0.2)
		if(missing(plots)) plots=TRUE
		
		if(missing(save)) save=NULL			
		
		HBSTM=hbstm.fit(HBSTM=HBSTM,nIter=nIter,nBurn=nBurn,time=control$time,timerem=control$timerem,plots=plots,posterior=posterior,save=save)
	}
	
	return(HBSTM)
	

	
}
setMethod(f="hbstm",signature=c("data.frame","matrix","data.frame","vector","vector","vector","ANY","ANY","ANY","ANY","logical","ANY","ANY","ANY","ANY"),definition=.HBSTM.creation)

.HBSTM.fit=function(HBSTM,nIter,nBurn,time,timerem,plots,posterior,save){
	
	if(missing(nIter)) nIter=1000
	if(missing(nBurn)) nBurn=as.integer(nIter*0.2)	
	if(missing(posterior)) posterior="mean"
	if(missing(time)) time=TRUE
	if(missing(timerem)) timerem=FALSE
	if(missing(plots)) plots=TRUE

	if(missing(save)){
		estYt.check(S=dim(HBSTM@newGrid)[1],T=dim(HBSTM@Zt)[2],nIter=nIter,nBurn=nBurn,val=mean(HBSTM@Zt))
		save=NULL
	}else{
		if(!is.null(save)){
			save.check(save=save,parameters=HBSTM@Parameters,nIter=nIter,nBurn=nBurn,S=dim(HBSTM@newGrid)[1],T=dim(HBSTM@Zt)[2],val=mean(HBSTM@Zt))
			HBSTM@MCMCclass=save
		}
	}

	if((nIter-nBurn)!=0){
		estYt=array(mean(HBSTM@Zt),dim=c(dim(HBSTM@newGrid)[1],dim(HBSTM@Zt)[2],nIter-nBurn))
	}
	msevect=rep(NA,nIter)
	
	if(HBSTM["seed"]!=-1) set.seed(HBSTM["seed"])
	
	S=dim(HBSTM["newGrid"])[1]
	cols=length(levels(as.factor(HBSTM["newGrid"][,1])))
	rows=length(levels(as.factor(HBSTM["newGrid"][,2])))
	T=dim(HBSTM["Zt"])[2]
	m=dim(HBSTM["Zt"])[1]
	P=as.matrix(cbind(rep(1,S),HBSTM["newGrid"]))
	identS=matrix(0,ncol=S,nrow=S)
	diag(identS)=1
	

	#Creación matrices auxiliares de Mt
	mtAux=CosSinMatrix.creation(S=S,T=T,seas=c(apply(as.matrix(1:length(HBSTM["Parameters"]["Mt"]["seas"])),1,function(el){return(HBSTM["Parameters"]["Mt"]["seas"][[el]]["w"])})))
	
	#Creación matrices espaciales auxiliares
	spatC=SpatC.creation(S=S,cols=cols,rows=rows,lengthSpatPar=HBSTM["Parameters"]["Mu"]["spatialMu"]["lags"],dirs=HBSTM["Parameters"]["Mu"]["spatialMu"]["dirs"])
	spatH=SpatH.creation(S=S,cols=cols,rows=rows,lengthSpatPar=HBSTM["Parameters"]["Xt"]["AR"][[1]]["subdiag"]["lags"],dirs=HBSTM["Parameters"]["Xt"]["AR"][[1]]["subdiag"]["dirs"])

	lagmodmax=max(HBSTM["Parameters"]["Xt"]["templags"])
	auxXtmat=matrix(0,nrow=S,ncol=lagmodmax)
	
	timeInter=matrix(NA,ncol=2,nrow=length(HBSTM["Parameters"]["Xt"]["templags"])+1)
	timeInter[,1]=c(max(HBSTM["Parameters"]["Xt"]["templags"])+1,(max(HBSTM["Parameters"]["Xt"]["templags"])+1)-HBSTM["Parameters"]["Xt"]["templags"])
	timeInter[,2]=c(T,T-HBSTM["Parameters"]["Xt"]["templags"])	
	
	auxTemplags=HBSTM["Parameters"]["Xt"]["templags"]
	auxSubdiag0=list()
	for(i in 1:length(auxTemplags)){
		auxSubdiag0[[i]]=HBSTM["Hyperpriors"]["Xt0"]["AR0"][[i]]["subdiag0"]
	}

	
	###### Update HBSTM["Parameters"]
	for(iter in 1:nIter){
		auxtim=proc.time()
		
		######## Update Yt
		invsigma2Y=1/HBSTM@Parameters@sigma2Y
		sigYt=chol2inv(chol(1/HBSTM@Parameters@sigma2E*t(HBSTM@K)%*%HBSTM@K+invsigma2Y*identS))
		HBSTM@Parameters@Yt=sigYt%*%t(1/HBSTM@Parameters@sigma2E*t(HBSTM@Zt)%*%HBSTM@K+invsigma2Y*t(matrix(rep(HBSTM@Parameters@Mu@muvect,T),nrow=S)+HBSTM@Parameters@Mt@Mt+HBSTM@Parameters@Xt@Xt))+t(chol(sigYt))%*%matrix(rnorm(S*T),ncol=T,nrow=S)
		rm(sigYt)

		######## Update Xt
		Xtaux=cbind(auxXtmat,HBSTM@Parameters@Xt@X0,HBSTM@Parameters@Xt@Xt,auxXtmat)
		auxsigma2=matrix(0,ncol=S,nrow=S)
		auxH=list()
		invsigma2N=1/HBSTM@Parameters@Xt@sigma2N
		for(i in 1:length(auxTemplags)){
			auxsigma2=auxsigma2+t(HBSTM@Parameters@Xt@AR[[i]]@H)%*%HBSTM@Parameters@Xt@AR[[i]]@H
			auxH[[i]]=HBSTM@Parameters@Xt@AR[[i]]@H
		}
		names(auxH)=HBSTM@Parameters@Xt@templags

		# X0
		sigma2X00inv=chol2inv(chol(HBSTM@Hyperpriors@Xt0@sigma2X00))
		sigX0=chol2inv(chol(invsigma2N*auxsigma2+sigma2X00inv))
		cholsigX0=t(chol(sigX0))
		auxX0=rep(0,S)
		if(length(auxTemplags)>1){
			for(i in 1:length(auxTemplags)){
				auxX0=auxX0+t(Xtaux[,i]-rowSums(apply(as.matrix(auxTemplags[-i]),1,function(el){
					return(auxH[[as.character(el)]]%*%Xtaux[,i+(lagmodmax+1)-el])
				}))
				)%*%auxH[[i]]
			}
		}else{
			auxX0=c(auxH[[1]]%*%Xtaux[,1+(lagmodmax+1)])
		}

		HBSTM@Parameters@Xt@X0=sigX0%*%t(invsigma2N*auxX0+t(HBSTM@Hyperpriors@Xt0@X00)%*%sigma2X00inv)+cholsigX0%*%as.matrix(rnorm(S))
		Xtaux[,lagmodmax+1]=HBSTM@Parameters@Xt@X0
		
		# Xt		
		sigXt=chol2inv(chol(invsigma2Y*identS+invsigma2N*(identS+auxsigma2)))
		cholsigXt=t(chol(sigXt))%*%matrix(rnorm(S*T),nrow=S,ncol=T)
		auxYtMuMt=invsigma2Y*(HBSTM@Parameters@Yt-matrix(HBSTM@Parameters@Mu@muvect,ncol=T,nrow=S)-HBSTM@Parameters@Mt@Mt)
		sigXt=t(sigXt)
		for(tim in (lagmodmax+2):(T+lagmodmax)){
			auxXt=rowSums(apply(as.matrix(auxTemplags),1,function(el){
				return(auxH[[as.character(el)]]%*%Xtaux[,tim-el])
			}))
			if(length(auxTemplags)>1){
				for(i in 1:length(auxTemplags)){
					auxlag=auxTemplags[i]
					auxXt=auxXt+(Xtaux[,tim+auxlag]-rowSums(apply(as.matrix(auxTemplags[-i]),1,function(el){
						return(auxH[[as.character(el)]]%*%Xtaux[,tim-el+auxlag])
					}))
					)%*%auxH[[i]]
				}			
			}else{
				auxXt=auxXt+(Xtaux[,tim+auxTemplags[1]])%*%auxH[[1]]
			}

			Xtaux[,tim]=(auxYtMuMt[,tim-lagmodmax]+invsigma2N*auxXt)%*%sigXt+cholsigXt[,tim-lagmodmax]			
		}
		HBSTM@Parameters@Xt@Xt=Xtaux[,(lagmodmax+1):(T+lagmodmax)]
		rm(sigX0,cholsigX0,auxX0,sigXt,auxXt,auxYtMuMt,Xtaux,cholsigXt,auxH)		
		
		########### Actualización avect y a0vect
		auxH=list(NULL)
		auxXt=list(NULL)
		auxXt[[1]]=HBSTM@Parameters@Xt@Xt[,timeInter[1,1]:timeInter[1,2]]
		auxmat=HBSTM@Parameters@Xt@Xt[,timeInter[1,1]:timeInter[1,2]]
		for(i in 1:length(auxTemplags)){
			auxH[[i]]=HBSTM@Parameters@Xt@AR[[i]]@H
			diag(auxH[[i]])=0
			auxXt[[i+1]]=HBSTM@Parameters@Xt@Xt[,timeInter[i+1,1]:timeInter[i+1,2]]
			auxmat=auxmat-HBSTM@Parameters@Xt@AR[[i]]@H%*%auxXt[[i+1]]
		}
	
		for(i in 1:length(auxTemplags)){
			invsigma2A=1/HBSTM@Parameters@Xt@AR[[i]]@sigma2A
			sigavect=chol2inv(chol(invsigma2N*diag(rowSums(auxXt[[i+1]]^2))+invsigma2A*(identS-HBSTM@Parameters@Xt@AR[[i]]@spatialA@Cmat)))			
			aux1=rowSums((auxmat+HBSTM@Parameters@Xt@AR[[i]]@H%*%auxXt[[i+1]]-auxH[[i]]%*%auxXt[[i+1]])*auxXt[[i+1]])			
			HBSTM@Parameters@Xt@AR[[i]]@avect=as.matrix(sigavect%*%t(invsigma2N*aux1+invsigma2A*t(HBSTM@Parameters@Xt@AR[[i]]@a0vect)%*%(identS-HBSTM@Parameters@Xt@AR[[i]]@spatialA@Cmat))+t(chol(sigavect))%*%as.matrix(rnorm(S)))
			
			siga0L=chol2inv(chol(invsigma2A*t(P)%*%(identS-HBSTM@Parameters@Xt@AR[[i]]@spatialA@Cmat)%*%P+diag(1/diag(HBSTM@Hyperpriors@Xt0@AR0[[i]]@siga0L0))))
			HBSTM@Parameters@Xt@AR[[i]]@a0L=as.matrix(siga0L%*%t(invsigma2A*t(HBSTM@Parameters@Xt@AR[[i]]@avect)%*%(identS-HBSTM@Parameters@Xt@AR[[i]]@spatialA@Cmat)%*%P+HBSTM@Hyperpriors@Xt0@AR0[[i]]@a0L0%*%diag(1/diag(HBSTM@Hyperpriors@Xt0@AR0[[i]]@siga0L0)))+t(chol(siga0L))%*%as.matrix(rnorm(3)))
			HBSTM@Parameters@Xt@AR[[i]]@a0vect=P%*%HBSTM@Parameters@Xt@AR[[i]]@a0L	
			
			diag(HBSTM@Parameters@Xt@AR[[i]]@H)=HBSTM@Parameters@Xt@AR[[i]]@avect
			
			for(j in HBSTM@Parameters@Xt@AR[[i]]@spatialA@dirs){
				HBSTM@Parameters@Xt@AR[[i]]@spatialA[j]=MHSpat(media=HBSTM@Parameters@Xt@AR[[i]]@avect,media0=HBSTM@Parameters@Xt@AR[[i]]@a0vect,invtau2=invsigma2A,alpha0=HBSTM@Hyperpriors@Xt0@AR0[[i]]@spatialA0[paste(j,"0",sep="")][1],sigalpha0=HBSTM@Hyperpriors@Xt0@AR0[[i]]@spatialA0[paste(j,"0",sep="")][2],namedir=j,spatParam=HBSTM@Parameters@Xt@AR[[i]]@spatialA,MatABPT=spatC)	
				HBSTM@Parameters@Xt@AR[[i]]@spatialA@Cmat=CmatCreation(spatParam=HBSTM@Parameters@Xt@AR[[i]]@spatialA,SpatC=spatC)
			}
		}
		rm(auxH,auxmat)
		
		
		auxtim2=proc.time()
		####### Update 	E,W,N,S,SE,NW,SW,NE


		auxSeq=1:(S*S)
		JmatAR=list(NULL)
		auxSubdiag=list()
		auxH=list()
		for(i in 1:length(auxTemplags)){
			JmatAR[[i]]=matrix(0,ncol=S,nrow=S)
			auxSubdiag[[i]]=HBSTM@Parameters@Xt@AR[[i]]@subdiag
			auxH[[i]]=HBSTM@Parameters@Xt@AR[[i]]@H
		}
		dirs=HBSTM@Parameters@Xt@AR[[1]]@subdiag@dirs
		
		for(k in 1:length(dirs)){
			dir0=paste(dirs[k],"0",sep="")
			auxSum=auxXt[[1]]
			
			Jb=lapply(as.matrix(rep(0,length(auxTemplags))),function(el){return(el)})
			for(j in 1:length(auxSubdiag[[i]][dirs[k]])){
				for(i in 1:length(auxTemplags)){
					JmatAR[[i]][spatH[dirs[k]][[j]]]=1		
					Jb[[i]]=Jb[[i]]+(auxSubdiag[[i]][dirs[k]][j]*JmatAR[[i]])
					JmatAR[[i]][spatH[dirs[k]][[j]]]=0
				}						
			}
			aux=list(NULL)
			for(i in 1:length(auxTemplags)){
				val=auxH[[i]][auxSeq[unlist(spatH[dirs[k]])]]
				auxH[[i]][auxSeq[unlist(spatH[dirs[k]])]]=0		
				aux[[i]]=(Jb[[i]]+auxH[[i]])%*%auxXt[[i+1]]	
				auxSum=auxSum-aux[[i]]
				auxH[[i]][auxSeq[unlist(spatH[dirs[k]])]]=val
			}
			for(i in 1:length(auxTemplags)){
				auxSum=auxSum+aux[[i]]
				for(j in 1:length(auxSubdiag[[i]][dirs[k]])){		
					Jb2=0
					for(l in c(1:length(auxSubdiag[[i]][dirs[k]]))[-j]){
						JmatAR[[i]][spatH[dirs[k]][[l]]]=1
						Jb2=Jb2+(auxSubdiag[[i]][dirs[k]][l]*JmatAR[[i]])
						JmatAR[[i]][spatH[dirs[k]][[l]]]=0				
					}	
					val=auxH[[i]][auxSeq[unlist(spatH[dirs[k]][-j])]]
					auxH[[i]][auxSeq[unlist(spatH[dirs[k]][-j])]]=0		
					JmatAR[[i]][spatH[dirs[k]][[j]]]=1
					aux2JmatAR=JmatAR[[i]]%*%auxXt[[i+1]]
					JmatAR[[i]][spatH[dirs[k]][[j]]]=0					
					invsigbvect=1/((invsigma2N)*sum((aux2JmatAR)^2)+1/auxSubdiag0[[i]][dir0][2])				
					auxSubdiag[[i]][dirs[k]][j]=rnorm(1,mean=invsigbvect*((invsigma2N)*sum((auxSum-(Jb2+auxH[[i]])%*%auxXt[[i+1]])*(aux2JmatAR))+auxSubdiag0[[i]][dir0][1]/auxSubdiag0[[i]][dir0][2]),sd=sqrt(invsigbvect))
					auxH[[i]][auxSeq[unlist(spatH[dirs[k]][-j])]]=val
					auxH[[i]][spatH[dirs[k]][[j]]]=auxSubdiag[[i]][dirs[k]][j]							
				}	
				auxSum=auxSum-aux[[i]]
			}
			
			rm(aux,auxSum,aux2JmatAR,invsigbvect,Jb2)
		}
		for(i in 1:length(auxTemplags)){
			HBSTM@Parameters@Xt@AR[[i]]@subdiag=auxSubdiag[[i]]
			HBSTM@Parameters@Xt@AR[[i]]@H=auxH[[i]]
		}
		
		######## Update Mu
		invsigma2Mu=1/HBSTM@Parameters@Mu@sigma2Mu
		sigmuvect=chol2inv(chol(invsigma2Mu*(identS-HBSTM@Parameters@Mu@spatialMu@Cmat)+T*invsigma2Y*identS))
		HBSTM@Parameters@Mu@muvect=as.matrix(sigmuvect%*%t(invsigma2Mu*t(HBSTM@Parameters@Mu@mu0vect)%*%(identS-HBSTM@Parameters@Mu@spatialMu@Cmat)+invsigma2Y*rowSums(HBSTM@Parameters@Yt-HBSTM@Parameters@Mt@Mt-HBSTM@Parameters@Xt@Xt))+t(chol(sigmuvect))%*%as.matrix(rnorm(S)))
		rm(sigmuvect)

		sigmu0L=chol2inv(chol(invsigma2Mu*t(P)%*%(identS-HBSTM@Parameters@Mu@spatialMu@Cmat)%*%P+diag(1/diag(HBSTM@Hyperpriors@Mu0@sigmu0L0))))
		HBSTM@Parameters@Mu@mu0L=as.matrix(sigmu0L%*%t(invsigma2Mu*t(HBSTM@Parameters@Mu@muvect)%*%(identS-HBSTM@Parameters@Mu@spatialMu@Cmat)%*%P+HBSTM@Hyperpriors@Mu0@mu0L0%*%diag(1/diag(HBSTM@Hyperpriors@Mu0@sigmu0L0)))+t(chol(sigmu0L))%*%as.matrix(rnorm(3)))	
		HBSTM@Parameters@Mu@mu0vect=P%*%HBSTM@Parameters@Mu@mu0L	
		rm(sigmu0L)
		for(i in HBSTM@Parameters@Mu@spatialMu@dirs){
			HBSTM@Parameters@Mu@spatialMu[i]=MHSpat(media=HBSTM@Parameters@Mu@muvect,media0=HBSTM@Parameters@Mu@mu0vect,invtau2=invsigma2Mu,alpha0=HBSTM@Hyperpriors@Mu0@spatialMu0[paste(i,"0",sep="")][1],sigalpha0=HBSTM@Hyperpriors@Mu0@spatialMu0[paste(i,"0",sep="")][2],namedir=i,spatParam=HBSTM@Parameters@Mu@spatialMu,MatABPT=spatC)	
			HBSTM@Parameters@Mu@spatialMu@Cmat=CmatCreation(spatParam=HBSTM@Parameters@Mu@spatialMu,SpatC=spatC)
		}
		
		######## Update Mt
		auxmat=HBSTM@Parameters@Yt-(matrix(HBSTM@Parameters@Mu@muvect,ncol=T,nrow=S)+HBSTM@Parameters@Xt@Xt)
		for(i in 1:length(HBSTM@Parameters@Mt@seas)){
			val=c(1:length(HBSTM@Parameters@Mt@seas))[-i]
			auxf0=c(HBSTM@Parameters@Mt@seas[[i]]@gvect)*mtAux@seas[[i]]@sinMat
			for(j in val){
				auxf0=auxf0+mtAux@seas[[j]]@cosMat*c(HBSTM@Parameters@Mt@seas[[j]]@fvect)+mtAux@seas[[j]]@sinMat*c(HBSTM@Parameters@Mt@seas[[j]]@gvect)
			}
			auxmeanf0=colSums(((t(auxmat-auxf0))%*%P)*mtAux@seas[[i]]@littlecosMat)	
			sigf0=chol2inv(chol(diag(1/diag(HBSTM@Hyperpriors@Mt0@seas0[[i]]@sigf0L0))+invsigma2Y*sum(cos(HBSTM@Parameters@Mt@seas[[i]]@w*c(1:T))^2)*t(P)%*%P))
			HBSTM@Parameters@Mt@seas[[i]]@f0L=sigf0%*%t(t(HBSTM@Hyperpriors@Mt0@seas0[[i]]@f0L0)%*%diag(1/diag(HBSTM@Hyperpriors@Mt0@seas0[[i]]@sigf0L0))+invsigma2Y*t(auxmeanf0))+t(chol(sigf0))%*%matrix(rnorm(3*1),ncol=1,nrow=3)
			HBSTM@Parameters@Mt@seas[[i]]@fvect=P%*%HBSTM@Parameters@Mt@seas[[i]]@f0L

			
			auxg0=auxf0-c(HBSTM@Parameters@Mt@seas[[i]]@gvect)*mtAux@seas[[i]]@sinMat+mtAux@seas[[i]]@cosMat*c(HBSTM@Parameters@Mt@seas[[i]]@fvect)
			auxmeang0=colSums(((t(auxmat-auxg0))%*%P)*mtAux@seas[[i]]@littlesinMat)
			sigg0=chol2inv(chol(diag(1/diag(HBSTM@Hyperpriors@Mt0@seas0[[i]]@sigg0L0))+invsigma2Y*sum(sin(HBSTM@Parameters@Mt@seas[[i]]@w*c(1:T))^2)*t(P)%*%P))
			HBSTM@Parameters@Mt@seas[[i]]@g0L=sigg0%*%t(t(HBSTM@Hyperpriors@Mt0@seas0[[i]]@g0L0)%*%diag(1/diag(HBSTM@Hyperpriors@Mt0@seas0[[i]]@sigg0L0))+invsigma2Y*t(auxmeang0))+t(chol(sigg0))%*%matrix(rnorm(3*1),ncol=1,nrow=3)
			HBSTM@Parameters@Mt@seas[[i]]@gvect=P%*%HBSTM@Parameters@Mt@seas[[i]]@g0L
			
			rm(val,auxf0,auxmeanf0,sigf0,auxg0,auxmeang0,sigg0)			
		}
		rm(auxmat)
		
		HBSTM@Parameters@Mt@Mt=matrix(0,ncol=T,nrow=S)
		for(i in 1:length(HBSTM@Parameters@Mt@seas)){
			HBSTM@Parameters@Mt@Mt=HBSTM@Parameters@Mt@Mt+matrix(HBSTM@Parameters@Mt@seas[[i]]@fvect,ncol=T,nrow=S)*mtAux@seas[[i]]@cosMat+
				matrix(HBSTM@Parameters@Mt@seas[[i]]@gvect,ncol=T,nrow=S)*mtAux@seas[[i]]@sinMat
		}
		
		


			
		######## Update sigma2s
		
		auxbetaN=auxXt[[1]]
		for(i in 1:length(auxTemplags)){
			HBSTM@Parameters@Xt@AR[[i]]@sigma2A=1/rgamma(1,S/2+HBSTM@Hyperpriors@Xt0@AR0[[i]]@sigma2A0[1],1/HBSTM@Hyperpriors@Xt0@AR0[[i]]@sigma2A0[2]+1/2*t(HBSTM@Parameters@Xt@AR[[i]]@avect-HBSTM@Parameters@Xt@AR[[i]]@a0vect)%*%(identS-HBSTM@Parameters@Xt@AR[[i]]@spatialA@Cmat)%*%(HBSTM@Parameters@Xt@AR[[i]]@avect-HBSTM@Parameters@Xt@AR[[i]]@a0vect))
			auxbetaN=auxbetaN-HBSTM@Parameters@Xt@AR[[i]]@H%*%auxXt[[i+1]]
		}
		HBSTM@Parameters@Xt@sigma2N=1/rgamma(1,S*T/2+HBSTM@Hyperpriors@Xt0@sigma2N0[1],1/HBSTM@Hyperpriors@Xt0@sigma2N0[2]+1/2*sum(auxbetaN^2))	
		rm(auxbetaN,auxXt)
		

		HBSTM@Parameters@Mu@sigma2Mu=1/rgamma(1,S/2+HBSTM@Hyperpriors@Mu0@sigma2Mu0[1],1/HBSTM@Hyperpriors@Mu0@sigma2Mu0[2]+1/2*t(HBSTM@Parameters@Mu@muvect-HBSTM@Parameters@Mu@mu0vect)%*%(identS-HBSTM@Parameters@Mu@spatialMu@Cmat)%*%(HBSTM@Parameters@Mu@muvect-HBSTM@Parameters@Mu@mu0vect))

		HBSTM@Parameters@sigma2Y=1/rgamma(1,S*T/2+HBSTM@Hyperpriors@sigma2Y0[1],1/HBSTM@Hyperpriors@sigma2Y0[2]+1/2*sum((HBSTM@Parameters@Yt-matrix(HBSTM@Parameters@Mu@muvect,nrow=S,ncol=T)-HBSTM@Parameters@Mt@Mt-HBSTM@Parameters@Xt@Xt)^2))
		
		HBSTM@Parameters@sigma2E=1/rgamma(1,m*T/2+HBSTM@Hyperpriors@sigma2E0[1],1/HBSTM@Hyperpriors@sigma2E0[2]+1/2*sum((HBSTM@Zt-HBSTM@K%*%HBSTM@Parameters@Yt)^2))


			
		if(iter>nBurn){
			estYt[,,iter-nBurn]=HBSTM@Parameters@Yt
			if(!is.null(save)){
				if(save=="Parameters"){
					HBSTM@MCMCsamp[[iter-nBurn]]=HBSTM@Parameters
				}else{
					Param=HBSTM@Parameters
					if(save=="Mu"){
						Param@Mt=new(Class="Mt")
						Param@Xt=new(Class="Xt")
						HBSTM@MCMCsamp[[iter-nBurn]]=Param
					}else{
						if(save=="Mt"){
							Param@Mu=new(Class="Mu")
							Param@Xt=new(Class="Xt")
							HBSTM@MCMCsamp[[iter-nBurn]]=Param
						}else{
							if(save=="Xt"){
								Param@Mu=new(Class="Mu")
								Param@Mt=new(Class="Mt")
								HBSTM@MCMCsamp[[iter-nBurn]]=Param
							}
						}
					}
				}
			}
		}

		if(time){
			if(iter==1){
				auxtim=proc.time()-auxtim
				cat("-------- The approximated execution time is: ", auxtim[3]*nIter/3600," hours------\n",sep="")
				
			}
		}
		aux=sum((HBSTM@Zt-(HBSTM@K%*%HBSTM@Parameters@Yt+HBSTM@Parameters@sigma2E))^2)
		cat("Iteration ",iter,": MSE = ",aux,"\n",sep="")
		msevect[iter]=aux
		if(timerem){
			if(iter>1){
				auxtim=proc.time()-auxtim
				cat("-------- The approximated time remaining is: ", auxtim[3]*(nIter-iter)/3600," hours------\n",sep="")
			}
		}

		if(plots){
			point=sample(1:dim(HBSTM@newGrid)[1],1)
			
			par(mfrow=c(2,2))
			plot(ts(msevect[1:iter]),main="MSE of the execution",ylab="MSE",xlab="iterations")
			aux=HBSTM@K%*%HBSTM@Parameters@Yt
			plot(ts(HBSTM@Zt[point,]),ylim=c(min(HBSTM@Zt[point,]),max(HBSTM@Zt[point,])+0.30*(max(HBSTM@Zt[point,])-min(HBSTM@Zt[point,]))),main=paste("Zt vs K*Yt in ",point," en iteration ",iter,sep=""),ylab="Values")
			lines(ts(aux[point,]),col=2)	
			legend("topright",col=c(1,2),lty=2,legend=c("Zt","K*Yt"))
			acf(HBSTM@Zt[point,]-aux[point,],lag.max=200,col=c(rep(1,8-1),2),main=paste("ACF of the residuals in point ",point,sep=""))	
			pacf(HBSTM@Zt[point,]-aux[point,],lag.max=200,col=c(rep(1,8-1),2),main=paste("PACF of the residuals in point ",point,sep=""))		
		}

		
	}

	if((nIter-nBurn)!=0){
		if(posterior=="mean"){
			HBSTM@fitted[,,1]=apply(estYt,c(1,2),mean)
			HBSTM@fitted[,,2]=apply(estYt,c(1,2),sd)
		}else{
			HBSTM@fitted[,,1]=apply(estYt,c(1,2),median)
			HBSTM@fitted[,,2]=apply(estYt,c(1,2),function(el){
				quantile(el,0.05)
			})	
			HBSTM@fitted[,,3]=apply(estYt,c(1,2),function(el){
				quantile(el,0.95)
			})			
		}
	}
	
	HBSTM@residuals=HBSTM@Zt-(HBSTM@K%*%HBSTM@fitted[,,1]+HBSTM@Parameters@sigma2E)
	HBSTM@iterations=HBSTM@iterations+nIter
	HBSTM@mse=c(HBSTM@mse,msevect)
	return(HBSTM)
}
setMethod(f="hbstm.fit",signature=c("HBSTM","ANY","ANY","ANY","ANY","ANY","ANY","ANY"),definition=.HBSTM.fit)





