.between.nptest <- function(data, perms=5000,  tail = NULL,statTest="t",testType="permutation",otherParams){
	
	if(is.null(statTest)) 
		if(is.null(otherParams$separatedX) || otherParams$separatedX)
			statTest="t" else
			statTest="F"
	
	if(!is.null(otherParams$rotationTest)) {
    if(otherParams$rotationTest) { testType="rotation"; rotationTest=TRUE} else rotationTest=FALSE
	} else 
	
	if(is.function(statTest)) {
		test <- statTest
	} else
	if(statTest=="t"){
		test <- .t.between.nptest		
	} else  
	if(statTest=="F"){ #ANOVAtype test, 1 column for each column of Y summarizing the dependence with all Xs			
	    test <- .F.between.nptest
	} else  
	if(statTest=="Trace"){ #ANOVAtype test, 1 column for each column of Y summarizing the dependence with all Xs			
		test <- .trace.between.nptest
	} else {warning("This test is not valid, nothing done."); test <- function() return()}

  environment(test) <- sys.frame(sys.nframe())
  out <- sys.frame(sys.nframe())
  return(out)
}

###########################################
.t.between.nptest <- function(){
  data=.getW(data)
  if((!is.null(otherParams$alsoMANOVA)&&otherParams$alsoMANOVA)|(!is.null(otherParams$onlyMANOVA)&&otherParams$onlyMANOVA)){  
    combData=.getLinCombAndWeights(data, otherParams)
    if(!is.null(otherParams$onlyMANOVA)&&otherParams$onlyMANOVA){
      data$Y =combData$Y
      data$W =combData$W
      nVar=combData$nVar
#       data$Su =combData$Su
    } else{
      data$Y =cbind(data$Y,comb=combData$Y)  
      data$W =cbind(data$W,comb=combData$W)
      nVar=c(rep(1,ncol(data$Y)),combData$nVar)
#       data$Su =rbind(cbind(data$Su,comb=matrix(NA,nrows(data$Su),ncol(combData$Su))),
#                      cbind(matrix(NA,nrows(data$Su),ncol(combData$Su)),comb=combData$Su))
    }
  }
	#colnames(data$W)=colnames(data$Y)
    perms <- make.permSpace(1:(nrow(data$Y) - max(0,ncol(data$Z))),perms,return.permIDs=TRUE,testType=testType)
	uni.test <- function(i,data){	
		data <- .orthoZ(list(X=data$W[,i]*data$X,Y=data$W[,i]*data$Y[,i,drop=FALSE],Z=data$W[,i]*data$Z,intercept = FALSE))
		#environment(.prod2t) <- sys.frame(sys.parent())
   	permT=.prod.perms(data,perms,testType=testType)
		permT=.prod2t(permT,data)
		permT
	}
  nX=ncol(data$X)
  permT = matrix(, perms$B,ncol(data$Y)*nX)
  for(i in 1:ncol(data$Y)) {
    permT[,(i-1)*nX+(1:nX)]=uni.test(i,data)
  }
  colnames(permT)=.getTNames(data$Y,if(nX>1) data$X else NULL)
  rownames(permT)=.getTRowNames(permT)
  if(!exists("nVar"))
    return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="w-t"))) else
      return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(nVar=nVar, Test="w-t")))
}
	
	
# #########################################
.trace.between.nptest <- function(){
	if(ncol(data$Y)==1) { warning("Only one response variable is inputed, use F statistic. Trace test is not performed.")
		return(NULL)
	} else {
	  nVar=NULL
	if(!is.null(otherParams$subsets)){
	  permT=NULL
    #si potrebbe tentare di fissare le dimensioni della matrice permT per velocizzare un po'
	    otherParams2=otherParams
	    otherParams2$subsets=NULL
	    for(i in 1:length(otherParams$subsets))
	    {
	      dataSub=.subsetData(data,subset=otherParams$subsets[[i]])
	      permT=cbind(permT,.trace.sim(dataSub,perms))
	      nVar=c(nVar,ncol(dataSub$Y))
	    }
	  colnames(permT)=names(otherParams$subsets)
	  } else {permT=.trace.sim(data,perms)
	          nVar=ncol(data$Y)}

	rownames(permT)=.getTRowNames(permT)
	return(list(permT=permT,perms=list(B=perms),tail=1,extraInfoPre=list(Test="F-trace",nVar=nVar)))
	}
}

.trace.sim <- function(data,perms){
  covs=data$covs
  for(i in 1:nrow(data$covs)) covs[i,,]=data$covs[i,,] +data$Su
  perms <- make.permSpace(covs,perms,testType="Simulation")	
  
  if(is.null(data$Z) || (ncol(data$Z)==0)) 
    {
    PZXPZ= data$X%*%solve(t(data$X)%*%data$X)%*%t(data$X)
    HZX= diag(nrow(PZXPZ)) - PZXPZ
    dfratio=(nrow(data$Y)- ncol(data$X))/ncol(data$X)
    }  else { 
       PZ= data$Z%*%solve(t(data$Z)%*%data$Z)%*%t(data$Z)
       PZX=cbind(data$Z,data$X)%*%solve(t(cbind(data$Z,data$X))%*%cbind(data$Z,data$X))%*%t(cbind(data$Z,data$X))
       PZXPZ= PZX-PZ
       HZX= diag(nrow(PZX)) - PZX
       rm(PZ,PZX)
       dfratio=(nrow(data$Y)- ncol(data$X) -ncol(data$Z) )/ncol(data$X)
     }
	.stat <- function(y) .tr(t(y)%*%PZXPZ %*% y)/.tr(t(y)%*%HZX%*%y) * dfratio
	permT=matrix(,perms$B,1)
  permT[1,]=.stat(data$Y)  
	for(i in 1:(perms$B-1)) permT[i+1,]=.stat(perms$rotFunct())
	colnames(permT)="F-trace"
	permT
}

#################################################
.F.between.nptest <- function(){
  data=.getW(data)
  if((!is.null(otherParams$alsoMANOVA)&&otherParams$alsoMANOVA)|(!is.null(otherParams$onlyMANOVA)&&otherParams$onlyMANOVA)){  
    combData=.getLinCombAndWeights(data, otherParams)
    if(!is.null(otherParams$onlyMANOVA)&&otherParams$onlyMANOVA){
      data$Y =combData$Y
      data$W =combData$W
      nVar=combData$nVar
      #       data$Su =combData$Su
    } else{
      data$Y =cbind(data$Y,comb=combData$Y)  
      data$W =cbind(data$W,comb=combData$W)
      nVar=c(rep(1,ncol(data$Y)),combData$nVar)
    }
  }
  
  perms <- make.permSpace(1:(nrow(data$Y) - max(0,ncol(data$Z))),perms,return.permIDs=TRUE,testType=testType)
  uni.test <- function(i,data){	
    data <- .orthoZ(list(X=data$W[,i]*data$X,Y=data$W[,i]*data$Y[,i,drop=FALSE],Z=data$W[,i]*data$Z,intercept = FALSE))
    P=.get.eigenv.proj.mat(data)
    permT=.prod.perms.P(data,perms,testType=testType,P)
    permT=.prod2F(permT,data)
    permT
  }
  permT = matrix(,perms$B,ncol(data$Y))
  for(i in 1:ncol(data$Y)) permT[,i]= uni.test(i,data)
  colnames(permT)=.getTNames(data$Y)
  rownames(permT)=.getTRowNames(permT)
  if(!exists("nVar"))
    return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="w-t"))) else
      return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(nVar=nVar, Test="w-t")))
}


#####################################
.getLinCombAndWeights<- function(data, otherParams){
  if(!is.null(otherParams$subsets)){
    dataComb=list(Y=NULL,W=NULL,nVar=NULL)
    otherParams2=otherParams
    otherParams2$subsets=NULL
    data=.getW(data)
    for(i in 1:length(otherParams$subsets))
      {
      dataSub=.subsetData(data,subset=otherParams$subsets[[i]])
      
      temp=.getLinCombAndWeights(dataSub,otherParams2)
      colnames(temp$Y)=paste(names(otherParams$subsets)[i],sep=":", colnames(temp$Y))
      colnames(temp$W)=paste(names(otherParams$subsets)[i],sep=":", colnames(temp$W))
      dataComb$Y=cbind(dataComb$Y,temp$Y)
      dataComb$W=cbind(dataComb$W,temp$W)
      dataComb$nVar=c(dataComb$nVar,temp$nVar)
    }
    return(dataComb)
  }
    
  if(ncol(data$Y)==1) 
    {W=as.matrix(1/sqrt(data$Su[1,1] + (data$covs))) 
     return(list(Y=data$Y,W=W,nVar=1))#Su=data$Su[1,1]))
  }
  
	if(is.null(otherParams$fastSumCombination) || (otherParams$fastSumCombination)){
		if(is.null(otherParams$linComb)) {
			if(is.null(otherParams$whichPCs)) whichPCs=1 else whichPCs=otherParams$whichPCs
      # standard deviations
      sd2s=diag(apply(data$covs,c(2,3),mean))+diag(data$Su)
			linComb=1/sqrt(sd2s)
			CORR=array(t(apply(data$covs,1,function(cov){ diag(linComb)%*%(cov+data$Su)%*%diag(linComb)})),dim(data$covs))
			CORR=apply(CORR,c(2,3),mean)
			linComb = (linComb)*prcomp(CORR)$rotation[,whichPCs,drop=FALSE]
      if(is.null(otherParams$onlyMANOVA)||!(otherParams$onlyMANOVA)) 
        linComb = cbind(linComb,sum=1)
      } else {
        if(length(otherParams$linComb)==1) otherParams$linComb=rep(otherParams$linComb,ncol(data$Y))
        linComb=otherParams$linComb
		}
		Y =data$Y%*%linComb
		W=array(,dim(Y))
		dimnames(W)=dimnames(Y)  
		for(j in 1:nrow(data$covs)) {
		  cov= data$Su + data$covs[j,,]
      W[j,]=1/sqrt(diag(t(linComb)%*%cov%*%linComb))
		}
		SuMulti=t(linComb)%*%data$Su%*%linComb
    } else{ #slow combination method
      if(is.null(otherParams$linComb)) {
        if(is.null(otherParams$whichPCs)) whichPCs=1:(1+(ncol(data$Y)==2)) else 
          whichPCs=otherParams$whichPCs
        meanCov=apply(data$covs,c(2,3),mean)
        linComb= diag(1/sqrt(diag(meanCov)))%*%prcomp(meanCov,scale. = TRUE)$rotation[,whichPCs,drop=FALSE]
        if(is.null(otherParams$onlyMANOVA)||!(otherParams$onlyMANOVA)) linComb = cbind(linComb,sum=1)
        } else {
          if(length(otherParams$linComb)==1) otherParams$linComb=rep(otherParams$linComb,ncol(data$Y))
          linComb=otherParams$linComb
        }
      Y=data$Y%*%linComb
      covMulti=array(t(apply(data$covs,1,function(cov){ diag(t(linComb)%*%cov%*%linComb)})),c(nrow(data$covs),ncol(linComb)))
      SuMulti=sapply(1:ncol(Y),function(i) .estimateSuMultiILS(Y=Y[,i,drop=FALSE],Z=as.matrix(cbind(data$X,data$Z)), S=array(covMulti[,i,drop=FALSE],c(nrow(Y),1,1))))
      W=array(,dim(Y))
      dimnames(W)=dimnames(Y)  
      for(j in 1:nrow(data$covs)) {
        cov= data$Su + data$covs[j,,]
        W[j,]=1/sqrt(covMulti[j,]+SuMulti)
      }	
	}
  if(is.null(dim(linComb)))
    nVar=length(linComb) else
      nVar=apply(linComb,2,function(x)sum(x!=0))
list(Y=Y,W=W,nVar=nVar)#,Su=SuMulti)
}


.subsetData <- function(data,subset){
  if(is.character(subset))  {
    sub=as.vector(unlist(sapply(subset,function(x) which(colnames(data$Y)==x))))
    if(length(sub)==0) warning(paste("no data found for subset with element: ",subset,sep="",". This will likely produce an error."))
    subset=sub
  }
  dataSub=data
  dataSub$Y=dataSub$Y[,subset,drop=FALSE]
  dataSub$Su=dataSub$Su[subset,subset,drop=FALSE]
  dataSub$covs=dataSub$covs[,subset,subset,drop=FALSE]
  dataSub$df.mod=dataSub$df.mod[,subset,drop=FALSE]
  dataSub
}
