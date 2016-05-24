flipMix <- function(modelWithin,X=NULL,Z=NULL,units, perms=1000, data=NULL, tail=NULL,
                    statTest=NULL,flipReturn, testType="permutation", 
                    Su=NULL, equal.se=FALSE,se=NA,replaceNA.coeffWithin="coeffMeans",
                    replaceNA.coeffWithin.se=replaceNA.coeffWithin, ...) {

  otherParams= list(...)
  if(is.null(otherParams$alsoMANOVA)) otherParams$alsoMANOVA=FALSE
  
  if(missing(flipReturn)||is.null(flipReturn)) 
  flipReturn=list(permT=TRUE,permP=FALSE,data=TRUE,call.env=TRUE)
  
  if(is.null(statTest) ) if(is.null(otherParams$separatedX)   || otherParams$separatedX)  
    { statTest="t" } else statTest="F"
  
  if(!(testType%in%c("permutation","rotation","simulation","symmetry"))) {
    if(is.null(otherParams$rotationTest) || (!otherParams$rotationTest) ) 
	{testType="permutation" } else { testType="rotation"}
  }
  
  if(is.null(statTest)) statTest="t"
  
  # store the call
  call <- match.call()
  if(!(is.list(data) && (!is.data.frame(data)))) {
     data<-obs2coeffWithin(modelWithin,X=X,Z=Z,units=units, data=data,equal.se=equal.se,se=se,
                        replaceNA.coeffWithin=replaceNA.coeffWithin,replaceNA.coeffWithin.se=replaceNA.coeffWithin.se,...)
 
  }
  rm(Z,X,modelWithin)
  #########
	N = nrow(data$coeffWithin)
	p = ncol(data$coeffWithin)
	############################## Estimate of random effects

	if(is.null(data$covs)) {data$covs=array(,c(N,p,p)); for( id in 1:N) data$covs[id,,]=diag(data$se[id,]^2)}
	if(is.null(Su)){
#     if(testType=='symmetry') #dalle simulazioni il guadagno in tempo pare essere minimo mentre il guadagno in potenza pare non irrilevante
#       {data$Su=sapply(1:ncol(data$coeffWithin), function(i)
#                                  .estimateSuMultiILS(Y=data$coeffWithin[,i,drop=FALSE],
#                                                      Z=as.matrix(cbind(data$X,data$Z)), 
#                                                      S=data$covs[,i,i,drop=FALSE])[1])
#        if(length(data$Su)>1) data$Su=diag(data$Su) else data$Su=matrix(data$Su)
#        }      else 
        data$Su=.estimateSuMultiILS(Y=data$coeffWithin,Z=as.matrix(cbind(data$X,data$Z)), S=data$covs)
	 } else {
		data$Su=Su; rm(Su) 
		}
	names(data)[names(data)=="coeffWithin"]="Y"
		if(length(unique(unlist(data$X)))>1){ # if X is not a constant perform dependence.nptest
		  res=.between.nptest(data, perms=perms, statTest=statTest[1], tail = tail, testType=testType,otherParams)
		  out=res$test()
		  #una pezza estetica:
		  if(ncol(data$X)==1)
		    colnames(out$permT)=paste(colnames(out$permT), "_|_",colnames(data$X),sep="")
		  res=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=res)
      if(length(statTest)>1){
        for(i in 2:length(statTest)) {
          ressub=.between.nptest(data, perms=perms, statTest=statTest[i], tail = tail, testType=testType,otherParams)
          out=ressub$test()
          #una pezza estetica:
          if(ncol(data$X)==1)
            colnames(out$permT)=paste(colnames(out$permT), "_|_",colnames(data$X),sep="")
          out=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=ressub)
          res=cFlip(res,out)          
        }
      }
		
      } else{ #otherwise perform a symmetry test
		  #estTypeWithin=c("none","H0","H1")
		  #type=H0 return a vector of estimates, 
      data$W=data$Y
		  data$W[,]=NA
		  for(j in 1:nrow(data$covs)){
		    data$W[j,]=1/sqrt(diag(data$Su) + diag(data$covs[j,,]))
		  }
      res=.symmetry.nptest(data, perms=perms, statTest=statTest[1],  tail = tail,testType=testType,...)
		  out=res$test()
      out$extraInfoPre=cbind(est.Su=diag(data$Su),out$extraInfoPre)
		  res=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=res)
      
      if(length(statTest)>1){
			for (i in 2:length(statTest)){
				data$W=matrix(,dim(data$Y))
        dimnames(data$W)=dimnames(data$Y)
		    for(j in 1:nrow(data$covs)){  
          data$W[j,]=1/sqrt(diag(data$Su) + diag(data$covs[j,,]))
		    }
			  ressub=.symmetry.nptest(data, perms=perms, statTest=statTest[i],  tail = tail,testType="t",...)
				out=ressub$test()
				out=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=res)
        res=cFlip(res,out)
			}
		}
	}
return(res)
}