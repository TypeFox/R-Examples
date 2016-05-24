

##################################################################################################################
####################################      Parameters creation     ################################################
##################################################################################################################



Parameter.creation=function(reglag,seas,spatlags,S,T){

	#Spatial param
	dirs=c("alpha","beta","phi","theta")[which(spatlags>0)]
	spatial=new(Class="SpatParam",Cmat=matrix(0,ncol=S,nrow=S),lags=spatlags[which(spatlags>0)],dirs=dirs)
	for(i in 1:length(dirs)){
		spatial[dirs[i]]=rep(0,spatial["lags"][i])
	}

	
	#Mu param
	mu=new(Class="Mu",muvect=as.matrix(rep(0,S)),mu0vect=as.matrix(rep(0,S)),mu0L=as.matrix(rep(0,3)),sigma2Mu=0,spatialMu=spatial)
	
	#Mt param
	aux=list(NULL)
	for(i in 1:length(seas)){	
		aux[[i]]=new(Class="Seas",w=seas[i],fvect=rep(0,S),f0L=as.matrix(rep(0,3)),gvect=rep(0,S),g0L=as.matrix(rep(0,3)))
	}	
	mt=new(Class="Mt",Mt=matrix(0,ncol=T,nrow=S),seas=aux)
	rm(aux)
	
	#Xt param
	auxDir=rep(spatlags,each=2)	
	dirs=c("east","west","north","south","southeast","northwest","southwest","northeast")[which(auxDir>0)]
	subdiag=new(Class="VectSubdiag",lags=auxDir[which(auxDir>0)],dirs=dirs)
	for(i in 1:length(dirs)){
		subdiag[dirs[i]]=rep(0,subdiag["lags"][i])
	}

	aux=list(NULL)
	for(i in 1:length(reglag)){	
		aux[[i]]=new(Class="Autoregressive",avect=as.matrix(rep(0,S)),a0vect=as.matrix(rep(0,S)),a0L=as.matrix(rep(0,3)),spatialA=spatial,sigma2A=0,H=matrix(0,ncol=S,nrow=S),subdiag=subdiag)
	}	
	names(aux)=c(reglag)
	
	xt=new(Class="Xt",Xt=matrix(0,ncol=T,nrow=S),X0=matrix(0,ncol=1,nrow=S),sigma2N=0,AR=aux,templags=reglag)
	
	return(new(Class="Parameters",Mu=mu,Mt=mt,Xt=xt,sigma2Y=0,sigma2E=0,Yt=matrix(0,ncol=T,nrow=S)))
	
}




##################################################################################################################################
##################################                Hyperpriors             ########################################################
##################################################################################################################################


HyperStandard.creation=function(reglag,Nseas,spatlags,S){
	
	#Spatial param
	nameSpat=c("alpha","beta","phi","theta")
	aux=list(NULL)
	for(i in 1:length(spatlags)){
		if(spatlags[i]>0){
			aux[[i]]=matrix(c(0,1000),nrow=2,ncol=spatlags[i])
			rownames(aux[[i]])=c("Mean","sigma2")
			colnames(aux[[i]])=paste(nameSpat[i],1:spatlags[i])
		}else{
			aux[[i]]=NULL
		}
	}
	spatial0=new(Class="SpatParam0",alpha0=aux[[1]],beta0=aux[[2]],phi0=aux[[3]],theta0=aux[[4]])
	rm(aux)
	
	#Mu param
	mu0=new(Class="Mu0",mu0L0=rep(0,3),sigmu0L0=diag(rep(1000,3)),sigma2Mu0=c(0.01,1000),spatialMu0=spatial0)
	
	#Mt param
	aux=list(NULL)
	for(i in 1:Nseas){	
		aux[[i]]=new(Class="Seas0",f0L0=rep(0,3),sigf0L0=diag(rep(1000,3)),g0L0=rep(0,3),sigg0L0=diag(rep(1000,3)))
	}	
	mt0=new(Class="Mt0",seas0=aux)
	rm(aux)

	#Xt param
	nameSpat=c("east","west","north","south","southeast","northwest","southwest","northeast")
	spatlags=rep(spatlags,2,each=2)
	aux=list(NULL)
	for(i in 1:length(spatlags)){
		if(spatlags[i]>0){
			aux[[i]]=matrix(c(0,1000),nrow=2,ncol=spatlags[i])
			rownames(aux[[i]])=c("Mean","sigma2")
			colnames(aux[[i]])=paste(nameSpat[i],1:spatlags[i])
		}else{
			aux[[i]]=NULL
		}
	}
	subdiag0=new(Class="VectSubdiag0",east0=aux[[1]],west0=aux[[1]],north0=aux[[2]],south0=aux[[2]],southeast0=aux[[3]],northwest0=aux[[3]],southwest0=aux[[4]],northeast0=aux[[4]])
	rm(aux)
	aux=list(NULL)
	for(i in 1:length(reglag)){	
		aux[[i]]=new(Class="Autoregressive0",a0L0=rep(0,3),siga0L0=diag(rep(1000,3)),sigma2A0=c(0.01,1000),spatialA0=spatial0,subdiag0=subdiag0)
	}	
	names(aux)=c(reglag)
	sigma2X00=diag(rep(0.01,S))
	xt0=new(Class="Xt0",X00=matrix(0,ncol=1,nrow=S),sigma2X00=sigma2X00,AR0=aux,sigma2N0=c(0.01,1000))
	
	hyper=new(Class="Hyperpriors",Mu0=mu0,Mt0=mt0,Xt0=xt0,sigma2Y0=c(0,1000),sigma2E0=c(0.01,1000))
	
}



x0L.check=function(x0L,name,pos=NULL,leng=3){
	if(!is.null(pos)){
		name=paste(pos,name)
	}
	if(length(x0L)!=leng){
		stop(paste("The length of the",name,"has to be",leng))
	}
	if(sum(is.na(x0L))>0){
		stop(paste("The are some NA's in the",name))
	}	
}

sigma2x0L.check=function(sigma2x0L,name,pos=NULL,leng=3){
	if(!is.null(pos)){
		name=paste(pos,name)
	}
	if((dim(sigma2x0L)[1]!=leng)|(dim(sigma2x0L)[2]!=leng)){
		stop(paste("The dimension of the",name,"has to be",leng,"X",leng))
	}
	if(sum(is.na(sigma2x0L))>0){
		stop(paste("The are some NA's in the",name))
	}	
}

sigma20.check=function(sigma2,name,pos=NULL){
	if(!is.null(pos)){
		name=paste(pos,name)
	}
	if(length(sigma2)!=3){
		stop(paste("The length of the",name,"has to be 2"))
	}
	if(sum(is.na(sigma2))>0){
		stop(paste("The are some NA's in the",name))
	}
}


spatial0.check=function(spatial,spatlags,name,pos=NULL){
	nameSpat=c("alpha","beta","phi","theta")
	if(!is.null(pos)){
		nameSpat=paste(pos,nameSpat)
	}
	nameSpat=paste(name,nameSpat)
	for(i in 1:length(spatlags)){
		if(spatlags[i]>0){
			if((dim(spatial[[i]])[1]!=2)|(dim(spatial[[i]])[2]!=spatlags[i])){
				stop(paste("The dimension of the",name,"has to be 2X",spatlags[i]))
			}	
			if(sum(is.na(spatial[[i]]))>0){
				stop(paste("The are some NA's in the",name))
			}			
		}else{
			if(!is.null(spatial[[i]])){
				stop(paste(nameSpat[[i]],"has to be NULL"))
			}
		}
	}	
}

subdiag0.check=function(subdiag0,spatlags,name,pos=NULL){
	nameSpat=c("east","west","north","south","southeast","northwest","southwest","northeast")
	if(!is.null(pos)){
		nameSpat=paste(pos,nameSpat)
	}
	nameSpat=paste(name,nameSpat)
	for(i in 1:length(spatlags)){
		if(spatlags[i]>0){
			if((dim(subdiag0[[i]])[1]!=2)|(dim(subdiag0[[i]])[2]!=spatlags[i])){
				stop(paste("The dimension of the",name,"has to be 2X",spatlags[i]))
			}	
			if(sum(is.na(subdiag0[[i]]))>0){
				stop(paste("The are some NA's in the",name))
			}			
		}else{
			if(!is.null(subdiag0[[i]])){
				stop(paste(nameSpat[[i]],"has to be NULL"))
			}
		}
	}	
}

Hyperpriors.control=function(hyperpriors,reglag,Nseas,spatlags,S){

	#Mu part
	x0L.check(x0L=hyperpriors["Mu0"]["mu0L0"],name="mu0L0")
	sigma2x0L.check(sigma2x0L=hyperpriors["Mu0"]["sigmu0L0"],name="sigmu0L0")
	sigma20.check(sigma2=hyperpriors["Mu0"]["sigma2Mu0"],name="sigma2Mu0")
	spatial0.check(spatial=hyperpriors["Mu0"]["spatialMu0"],spatlags=spatlags,name="Mu")

	#Mt part
	for(i in 1:Nseas){	
		x0L.check(x0L=hyperpriors["Mt0"]["Seas0"][[i]]["f0L0"],name=paste("f0L0",i,sep=""),pos=paste("seasonal",i))
		sigma2x0L.check(sigma2x0L=hyperpriors["Mt0"]["Seas0"][[i]]["sigf0L0"],name=paste("sigf0L0",i,sep=""),pos=paste("seasonal",i))
		x0L.check(x0L=hyperpriors["Mt0"]["Seas0"][[i]]["g0L0"],name=paste("g0L0",i,sep=""),pos=paste("seasonal",i))
		sigma2x0L.check(sigma2x0L=hyperpriors["Mt0"]["Seas0"][[i]]["sigg0L0"],name=paste("sigg0L0",i,sep=""),pos=paste("seasonal",i))
	}	

	#Xt part
	for(i in 1:length(reglag)){	
		x0L.check(x0L=hyperpriors["Xt0"]["AR0"][[i]]["a0L0"],name=paste("a0L0",i,sep=""),pos=paste("autoregressive",i))
		sigma2x0L.check(sigma2x0L=hyperpriors["Xt0"]["AR0"][[i]]["a0L0"],name=paste("siga0L0",i,sep=""),pos=paste("autoregressive",i))
		sigma20.check(sigma2=hyperpriors["Xt0"]["AR0"][[i]]["sigma2A0"],name="sigma2A0",pos=paste("autoregressive",i))
		spatial0.check(spatial=hyperpriors["Xt0"]["AR0"][[i]]["spatialA0"],spatlags=spatlags,name="Xt",pos=paste("autoregressive",i))
		subdiag0.check(subdiag0=hyperpriors["Xt0"]["AR0"][[i]]["subdiag0"],spatlags=rep(spatlags,2,each=2),name="Xt",pos=paste("autoregressive",i))
	}	
	x0L.check(x0L=hyperpriors["Xt0"]["X00"],name="X00",leng=S)
	sigma2x0L.check(sigma2x0L=hyperpriors["Xt0"]["sigma2X00"],name="sigma2X00",leng=S)
	sigma20.check(sigma2=hyperpriors["Xt0"]["sigma2N0"],name="sigma2N0")

	sigma20.check(sigma2=hyperpriors["sigma2Y0"],name="sigma2Y0")
	sigma20.check(sigma2=hyperpriors["sigma2E0"],name="sigma2E0")
	
	return(hyperpriors)
}



##################################################################################################################################
##################################               Initial values            #######################################################
##################################################################################################################################




initialStand.creation=function(Parameters,Zt,P){
	
	S=dim(P)[1]
	T=dim(Zt)[2]
	identS=diag(rep(1,S))
	#Mu part
	meanZt=apply(Zt,1,mean)
	meanMod=lm(meanZt~P[,2:3])
	coe=coef(meanMod)
	coe[is.na(coe)]=0
	Parameters["Mu"]["mu0L"]=as.matrix(coe)
	Parameters["Mu"]["sigma2Mu"]=sd(meanZt)^2
	Parameters["Mu"]["mu0vect"]=P%*%Parameters["Mu"]["mu0L"]
	invsigmuvect=solve(identS-Parameters["Mu"]["spatialMu"]["Cmat"])
	Parameters["Mu"]["muvect"]=as.matrix(mvrnorm(1,mu=Parameters["Mu"]["mu0vect"],Sigma=Parameters["Mu"]["sigma2Mu"]*invsigmuvect))
	rm(invsigmuvect)

	
	#Mt param
	meanZt=apply(Zt,2,mean)
	seasMat=NULL
	tim=1:dim(Zt)[2]
	meanP=matrix(apply(P,2,mean),ncol=3,nrow=dim(Zt)[2])
	for(i in 1:length(Parameters["Mt"]["seas"])){
		cosmat=meanP*matrix(cos(Parameters["Mt"]["seas"][[i]]["w"]*tim),ncol=3,nrow=dim(Zt)[2])
		sinmat=meanP*matrix(sin(Parameters["Mt"]["seas"][[i]]["w"]*tim),ncol=3,nrow=dim(Zt)[2])
		seasMat=cbind(seasMat,cosmat,sinmat)
	}
	
	seasMod=lm(meanZt~seasMat)
	coe=coef(seasMod)[-1]
	coe[is.na(coe)]=0
	for(i in 1:length(Parameters["Mt"]["seas"])){
		Parameters["Mt"]["seas"][[i]]["f0L"]=as.matrix(coe[1:3])
		Parameters["Mt"]["seas"][[i]]["g0L"]=as.matrix(coe[4:6])
		coe=coe[-c(1:6)]
		Parameters["Mt"]["seas"][[i]]["fvect"]=P%*%Parameters["Mt"]["seas"][[i]]["f0L"]
		Parameters["Mt"]["seas"][[i]]["gvect"]=P%*%Parameters["Mt"]["seas"][[i]]["f0L"]
		Parameters["Mt"]["Mt"]=Parameters["Mt"]["Mt"]+matrix(Parameters["Mt"]["seas"][[i]]["fvect"],ncol=T,nrow=S)*matrix(cos((2*pi/Parameters["Mt"]["seas"][[i]]["w"])*c(1:T)),ncol=T,nrow=S,byrow=TRUE)+matrix(Parameters["Mt"]["seas"][[i]]["gvect"],ncol=T,nrow=S)*matrix(sin((2*pi/Parameters["Mt"]["seas"][[i]]["w"])*c(1:T)),ncol=T,nrow=S,byrow=TRUE)	
	}	

	#Xt param
	Parameters["Xt"]["sigma2N"]=1
	for(i in 1:length(Parameters["Xt"]["templags"])){
		Parameters["Xt"]["AR"][[i]]["sigma2A"]=0.01
	}	
	
	#General
	Parameters["sigma2E"]=10
	Parameters["sigma2Y"]=10	
	return(Parameters)
}

difvectLeng=function(vect1,vect2,name){
	if(length(vect1)!=length(vect2)) stop(paste("The length of the initial value",name,"is different of the parameter"))
	return(NULL)
}

difmatDim=function(mat1,mat2,name){
	if((dim(mat1)[1]==dim(mat2)[1])&(dim(mat1)[2]==dim(mat2)[2])) stop(paste("The length of the initial value",name,"is different of the parameter"))
	return(NULL)
}

spatial.check=function(param,initial,spatdir,name){
	difvectLeng(param["dirs"],initial["dirs"],name=paste(name,"directions"))
	if(sum(param["dirs"]==initial["dirs"])!=length(param["dirs"])){
		stop(paste("The spatial directions of the initial value",name,"are different of the parameter"))
	}else{
		for(i in param["dirs"]){
			difvectLeng(param[i],initial[i],name=paste(name,"i"))
		}		
	}
	difmatDim(param["Cmat"],initial["Cmat"],name=paste(name,"Cmat"))
	return(NULL)
}


initial.control=function(Parameters,initialvalues){
	
	#Mu param
	difmatDim(Parameters["Mu"]["mu0L"],initialvalues["Mu"]["mu0L"],name="mu0L")
	difmatDim(Parameters["Mu"]["mu0vect"],initialvalues["Mu"]["mu0vect"],name="mu0vect")
	difmatDim(Parameters["Mu"]["muvect"],initialvalues["Mu"]["muvect"],name="muvect")
	difvectLeng(Parameters["Mu"]["sigma2Mu"],initialvalues["Mu"]["sigma2Mu"],name="mu0L")
	spatial.check(param=Parameters["Mu"]["spatialMu"],initial=initialvalues["Mu"]["spatialMu"],name="Mu")
	
	Parameters["Mu"]["mu0L"]=initialvalues["Mu"]["mu0L"]
	Parameters["Mu"]["mu0vect"]=initialvalues["Mu"]["mu0vect"]
	Parameters["Mu"]["muvect"]=initialvalues["Mu"]["muvect"]
	Parameters["Mu"]["sigma2Mu"]=initialvalues["Mu"]["sigma2Mu"]
	
	#Mt param
	
	for(i in 1:length(initialvalues["Mt"]["seas"])){
		difmatDim(Parameters["Mt"][[i]]["f0L"],initialvalues["Mt"][[i]]["f0L"],name=paste("seasonal",i,"f0L"))
		difmatDim(Parameters["Mt"][[i]]["g0L"],initialvalues["Mt"][[i]]["g0L"],name=paste("seasonal",i,"g0L"))
		difvectLeng(Parameters["Mt"][[i]]["fvect"],initialvalues["Mt"][[i]]["fvect"],name=paste("seasonal",i,"fvect"))
		difvectLeng(Parameters["Mt"][[i]]["gvect"],initialvalues["Mt"][[i]]["gvect"],name=paste("seasonal",i,"gvect"))
		
		Parameters["Mt"][[i]]["f0L"]=initialvalues["Mt"][[i]]["f0L"]
		Parameters["Mt"][[i]]["g0L"]=initialvalues["Mt"][[i]]["g0L"]
		Parameters["Mt"][[i]]["fvect"]=initialvalues["Mt"][[i]]["fvect"]
		Parameters["Mt"][[i]]["gvect"]=initialvalues["Mt"][[i]]["gvect"]
	}
	
	difmatDim(Parameters["Mt"]["Mt"],initialvalues["Mt"]["Mt"],name="Mt")
	Parameters["Mt"]=initialvalues["Mt"]
	
	#Xt param
	difmatDim(Parameters["Xt"]["Xt"],initialvalues["Xt"]["Xt"],name="Xt")
	difmatDim(Parameters["Xt"]["X0"],initialvalues["Xt"]["X0"],name="X0")
	difvectLeng(Parameters["Xt"]["sigma2N"],initialvalues["Xt"]["sigma2N"],name="sigma2N")
	
	Parameters["Xt"]["Xt"]=initialvalues["Xt"]["Xt"]
	Parameters["Xt"]["X0"]=initialvalues["Xt"]["X0"]
	Parameters["Xt"]["sigma2N"]=initialvalues["Xt"]["sigma2N"]
	
	for(i in 1:length(initialvalues["Xt"]["templags"])){
		
		difmatDim(Parameters["Xt"][[i]]["a0L"],initialvalues["Xt"][[i]]["a0L"],name=paste("autorregresive",i,"a0L"))
		difmatDim(Parameters["Xt"][[i]]["a0vect"],initialvalues["Xt"][[i]]["a0vect"],name=paste("autorregresive",i,"a0vect"))
		difmatDim(Parameters["Xt"][[i]]["avect"],initialvalues["Xt"][[i]]["avect"],name=paste("autorregresive",i,"avect"))
		difvectLeng(Parameters["Xt"][[i]]["sigma2A"],initialvalues["Xt"][[i]]["sigma2A"],name=paste("autorregresive",i,"sigma2A"))
		difmatDim(Parameters["Xt"][[i]]["H"],initialvalues["Xt"][[i]]["H"],name=paste("autorregresive",i,"H"))
		spatial.check(param=Parameters["Xt"][[i]]["spatialA"],initial=initialvalues["Xt"][[i]]["spatialA"],name=paste("autorregresive",i))
		spatial.check(param=Parameters["Xt"][[i]]["subdiag"],initial=initialvalues["Xt"][[i]]["subdiag"],name=paste("autorregresive",i))
	
		Parameters["Xt"][[i]]["a0L"]=initialvalues["Xt"][[i]]["a0L"]
		Parameters["Xt"][[i]]["a0vect"]=initialvalues["Xt"][[i]]["a0vect"]
		Parameters["Xt"][[i]]["avect"]=initialvalues["Xt"][[i]]["avect"]
		Parameters["Xt"][[i]]["sigma2A"]=initialvalues["Xt"][[i]]["sigma2A"]
		Parameters["Xt"][[i]]["H"]=initialvalues["Xt"][[i]]["H"]
		Parameters["Xt"][[i]]["spatialA"]=initial=initialvalues["Xt"][[i]]["spatialA"]
		Parameters["Xt"][[i]]["subdiag"]=initial=initialvalues["Xt"][[i]]["subdiag"]
		
	}	
	
	
	#General param
	difvectLeng(Parameters["sigma2Y"],initialvalues["sigma2Y"],name="sigma2Y")
	difvectLeng(Parameters["sigma2E"],initialvalues["sigma2E"],name="sigma2E")
	difmatDim(Parameters["Yt"],initialvalues["Yt"],name="Yt")

	Parameters["sigma2Y"]=initialvalues["sigma2Y"]
	Parameters["sigma2E"]=initialvalues["sigma2E"]
	Parameters["Yt"]=initialvalues["Yt"]
		
	return(Parameters)		
}



