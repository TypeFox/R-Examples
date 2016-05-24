obs2coeffWithin <- function(modelWithin,units,X=NULL, Z=NULL, data=NULL,equal.se=FALSE,se=NA,
                            replaceNA.coeffWithin=NA,replaceNA.coeffWithin.se=Inf,...){
  
  if(is(modelWithin,"formula"))
    modelWithin <- lm(modelWithin,data=data)
  
  if( is(units,"formula") && (!is.null(data)) ){
    units <- model.frame(units, data=data, drop.unused.levels = TRUE)
  }
  units=unlist(units)
  
  
  .obs2coeffWithin<-function(modelWithin,units=units, data=data,se=se,...){
		if(is(modelWithin,"vgam")|is(modelWithin,"vglm")) #e quindi anche vglm ed anche 
			out=.obs2coeffWithin.vglm(modelWithin,units=units, data=data,se=se,...) else
		if(is(modelWithin,"glm")) 
			out=.obs2coeffWithin.glm(modelWithin,units=units, data=data,se=se,...) else
		if(is(modelWithin,"lm")) #occhio che gli lm univariati sono anche glm (anche) in R
			out=.obs2coeffWithin.lm(modelWithin,units=units, data=data,se=se,...) 
		out
	}

	if(is(modelWithin,"list")){
		.out=list()
		for(i in 1:length(modelWithin)){
			.out[[i]]=.obs2coeffWithin(modelWithin[[i]],units=units, data=data,se=se,...) 
		}
		names(.out)=names(modelWithin)
		out=.unlist.obs2coeffWithin(.out)
	} else	out=.obs2coeffWithin(modelWithin,units=units, data=data,se=se,...) 
	
	########### missing values are replaced with average among subjects and other things to make up the output
	out=.fixDataREff(out,units=units,X=X,Z=Z,data=data,replaceNA.coeffWithin=replaceNA.coeffWithin,replaceNA.coeffWithin.se=replaceNA.coeffWithin.se,
                   equal.se=equal.se,se=se,ncoef=ncol(out$coeffWithin),...)
out
}



######################################
############# make Y (light version, version 1) for lm with possibly multivariate response
##############################
.obs2coeffWithin.lm <- function(modelWithin, units, data, se=NA, replaceNA.coeffWithin=NA,replaceNA.coeffWithin.se=NA,...){

	if (!is(modelWithin,"lm")) print("only lm allowed in .obs2coeffWithin.lm")
	
	idClust=unique(units)
	names(idClust)=idClust
	nUnit=length(unique(idClust))
	ncoef=length(coefficients(modelWithin))
  colNames= colnames(vcov(modelWithin))
# 	colNames =  paste(rep(colnames(coefficients(modelWithin)),
#                     rep(length(rownames(coefficients(modelWithin))),
#                         length(colnames(coefficients(modelWithin))))),
#                     rep(rownames(coefficients(modelWithin)),length(colnames(coefficients(modelWithin)))),sep=".")
  if(length(colNames)<1) colNames=names(coefficients(modelWithin))
	
	res=.getEmptyFlipMix(colNames = colNames,idClust=idClust,nUnit=nUnit,ncoef=ncoef)
	
	for(idsel in idClust )  {
  	newmodel = try(update(modelWithin,contrasts=modelWithin$contrasts,subset=(units==idsel),na.action=na.omit),silent = TRUE)
		if(!is(newmodel,"try-error")){
			res$coeffWithin[idsel,]=as.vector(coef(newmodel))
  	res$se[idsel,]=if(!is.null(dim(coefficients(modelWithin)))) 
			                   as.vector(sapply(summary(newmodel),
							      function(x) {out=rep(NA,length(rownames(modelWithin$coefficients))); 
								               names(out)=rownames(modelWithin$coefficients); 
											   out[names(x$coef[,"Std. Error"])]=x$coef[,"Std. Error"];
											   out})) else 
								as.vector(summary(newmodel)$coef[,"Std. Error"])
								
# 			temp=if(!is.null(dim(coefficients(modelWithin)))) as.vector(sapply(summary(newmodel),function(x)x$df[1:2])) else as.vector(summary(newmodel)$df[1:2])
			res$df.mod[idsel,] = 1
			res$df.res[idsel,] =  newmodel$df.residual
        temp=vcov(newmodel)
			res$covs[idsel,colnames(temp),colnames(temp)]= temp
#       res$se[idsel,colnames(temp)]=sqrt(diag(temp))
		}
	}
	####? serve?
	if(ncoef==1) { dimnames(res$covs)[[2]] <- dimnames(res$covs)[[3]] <- colNames}
	
	res
}


#############questa versione mette anche i nomi alle righe delle Y, se ecc
#############inoltre restituisce matrici di covarianze dei coefficients e non solo la diagonale
.obs2coeffWithin.glm <- function(modelWithin, units, se=NA, replaceNA.coeffWithin=NA,replaceNA.coeffWithin.se=NA,...){
### MODEL NON ? UNA LISTA

	idClust=unique(units)
	names(idClust)=idClust
	nUnit=length(unique(idClust))
	ncoef=length( coefficients(modelWithin)) 
	colNames = names(coefficients(modelWithin))

	res=.getEmptyFlipMix(colNames=colNames,idClust=idClust,ncoef=ncoef,nUnit=nUnit)
	
	for(idsel in idClust )  {
		newmodel = try(update(modelWithin,subset=(units==idsel)),silent = TRUE)
		if(is(newmodel,"try-error")){
				# temp<-matrix(NA,nrow=length(modelWithin$coef),ncol=2)
				# colnames(temp)<-c("Estimate","Std. Error")
		} else{	
				temp <- summary(newmodel)$coef
			
				#TODO renderla una matrice triangolare per risparmiare memoria? vale la pena?
				res$coeffWithin[idsel,rownames(temp)]=temp[,"Estimate"]
				res$se[idsel,rownames(temp)]=temp[,"Std. Error"]	()
				res$df.mod[idsel,rownames(temp)] = rep(dim(temp)[1],dim(temp)[1])
				res$df.res[idsel,1] =  sum(units==idsel)-dim(temp)[1]
				res$covs[idsel,rownames(temp),rownames(temp)]=summary(newmodel)$cov.scaled
				res$dispersion[idsel,1]=summary(newmodel)$dispersion
		}
	}
	
	res
}



#############questa versione mette anche i nomi alle righe delle Y, se ecc
#############inoltre restituisce matrici di covarianze dei coefficients e non solo la diagonale
.obs2coeffWithin.vglm <- function(modelWithin, units, se=NA, replaceNA.coeffWithin=NA,replaceNA.coeffWithin.se=NA,...){
#require(VGAM)
### MODEL NON ? UNA LISTA
	update.vglm <- function(modelWithin,subset){
		TEMP_DATA_REDUCTED=eval(modelWithin@call$data)[subset,]
		modelWithin@call$data=as.name("TEMP_DATA_REDUCTED")
		res=eval(modelWithin@call)
		res
	}
	

	idClust=unique(units)
	names(idClust)=idClust
	nUnit=length(unique(idClust))
	ncoef=length( coefficients(modelWithin))
	colNames=names(coefficients(modelWithin))
	res=.getEmptyFlipMix(colNames=colNames,idClust=idClust,ncoef=ncoef,nUnit=nUnit)
	
	for(idsel in idClust )  {
		newmodel = try(update.vglm(modelWithin,subset=(units==idsel)),silent = TRUE)
		if(is(newmodel,"try-error")){
				# temp<-matrix(NA,nrow=length(coefficients(modelWithin)),ncol=2)
				# colnames(temp)<-c("Estimate","Std. Error")
				# covCoeff=matrix(NA,nrow=dim(temp)[1],ncol=dim(temp)[1]) 
		} else{	
			temp <- summary(newmodel)@coef3
			colnames(temp)[colnames(temp)=="Value"]="Estimate"
			covCoeff <- summary(newmodel)@ cov.unscaled
			
				#TODO renderla una matrice triangolare per risparmiare memoria? vale la pena?
				res$coeffWithin[idsel,rownames(temp)]=temp[,"Estimate"]
				res$se[idsel,rownames(temp)]=temp[,"Std. Error"]											
				res$df.mod[idsel] = rep(nrow(temp),nrow(temp))
				res$df.res[idsel,rownames(temp)] =  sum(units==idsel)-dim(temp)[1]
				res$covs[idsel,rownames(temp),rownames(temp)]=covCoeff
		}
	}
	
	res
}




#######################################
.fixDataREff <- function(res,units,X=X,Z=Z,data,equal.se,se=NA,replaceNA.coeffWithin,replaceNA.coeffWithin.se,ncoef,...){
#replaceNA.coeffWithin.se e replaceNA.coeffWithin accettano "varMean", "colMean" o un real 
#per replaceNA.coeffWithin=NA e replaceNA.coeffWithin.se=NA, Y vale la media e se max(se)*1E3 per colonna  - equivale a missing coefficient
  if( is(units,"formula") && (!is.null(data)) )
    units <- model.frame(units, data, drop.unused.levels = TRUE)      
  N=length(unique(units))
  
  if( (is(X,"formula")||is(Z,"formula")) && (!is.null(data)) ){
    firstIDs <- sapply(unique(unlist(units)),function(x) which(unlist(units)==x)[1] )
    data <- data[firstIDs,,drop=FALSE]
    rownames(data)=rownames(res$coeffWithin)
  }
  
  if(!is.null(Z)){
    if(is(Z,"formula")) Z <-model.matrix(Z, data=data)  
    Z <- .getNull(Z, data=data, n=nrow(data))$Z
  }
  
  if(!is.null(X)){
    if(is(X,"formula")) X <-model.matrix(X, data=data)
    X <- .getAlternative(X, data=data, n=nrow(data),dummyfy=list(...)$dummyfy)
    if(!is.null(Z) && ncol(Z)>0)  X <- X[,!.getIntercept(X),drop=FALSE]
  } else
    X <- matrix(1,nrow=N)
  
  namesBetween=setdiff(unique(c(colnames(X), colnames(Z))),"(Intercept)")
  
  if(nrow(X)==length(as.vector(units))){
    XX=c()
    for(i in 1:ncol(X)){
      XX=cbind(XX, as.matrix(table(units,X[,i]))[,-1,drop=FALSE])
    }
    X=XX
    rm(XX)
    namesBetween=c(namesBetween,colnames(X))
  }
  
  if((!is.null(Z)) && (nrow(Z)==length(as.vector(units)))){
    ZZ=c()
    for(i in 1:ncol(Z))  {
      ZZ=cbind(ZZ,as.matrix(table(units,Z[,i]))[,-1,drop=FALSE])
    }
    Z=ZZ
    rm(ZZ)
    namesBetween=c(namesBetween,colnames(Z))
  }
  
  res$X=X  
  res$Z=Z
  rm(X,Z)
  
  res=.testedWithin(res,includeOnly= 
                      if(!is.null(list(...)$test.coeffWithin)) list(...)$test.coeffWithin else NULL
                      ,exclude=namesBetween)


	
  #replace estimated res$se with imposed (not missing) se. 
	if(!is.null(se))if(!is.na(se)){
	#if se has not the same size of Y build it as matrix
		if(!all(dim(as.matrix(se))==c(length(idClust),ncoef))) se=matrix(se[1:min(ncoef,prod(dim(se)))],ncol=ncoef,nrow=length(idClust),byrow=TRUE)
	res$se[!is.na(se)] <- se[!is.na(se)]
	}
	##NULL df.res are forces to 0
	if(!is.null(res$df.res)) res$df.res[is.na(res$df.res)]=0

	if(any(is.na(res$coeffWithin))){
	  if(!is.na(replaceNA.coeffWithin)){
      if(is.numeric(replaceNA.coeffWithin)){
        warning(paste("Some estimated coefficient within unit is NA, the value ",replaceNA.coeffWithin," will be imputed.",sep=""))
        res$coeffWithin[is.na(res$coeffWithin)]=rep(as.double(replaceNA.coeffWithin),
                                                    length.out=sum(is.na(res$coeffWithin)))
        
        } else {
          warning("Some estimated coefficient within unit is NA, the value will be imputed (estimated from the data).")
          if(replaceNA.coeffWithin=="coeffMeans"){
            for(i in which(apply(res$coeffWithin,2,function(x) any(is.na(x))) )) 
              res$coeffWithin[is.na(res$coeffWithin[,i]),i]=mean(res$coeffWithin[,i],na.rm=TRUE) 
            } else if(replaceNA.coeffWithin=="unitMeans") {
              for(i in which(apply(res$coeffWithin,1,function(x) any(is.na(x))) )) 
                res$coeffWithin[i,is.na(res$coeffWithin[i,])]=mean(res$coeffWithin[i,],na.rm=TRUE)
            }
        }
	  }
	}
	
	if((!is.null(equal.se)) && equal.se) {
		w=res$df.res
    if(ncol(w)!=ncol(res$se)) w=matrix(w,ncol=ncol(res$se),nrow=nrow(res$se))
		if(any(is.na(res$se))) {
			idNA=is.na(res$se)
			res$se[idNA]=0
			}
    if(is.null(res$dispersion)|| any(is.na(res$dispersion))){
      res$se= matrix(sqrt(apply(res$se^2 * w,2,sum)/apply(w,2,sum)),
                     byrow=TRUE,nrow=nrow(res$se),ncol=ncol(res$se))/sqrt(w)
      } else  { 
        if(ncol(res$dispersion)!=ncol(res$se)) 
          dispersion=matrix(res$dispersion,ncol=ncol(res$se),nrow=nrow(res$se)) else
            dispersion=res$dispersion
       res$se= matrix(sqrt(apply(res$se^2 *dispersion^2* w,2,sum)/apply(w,2,sum)),
                      byrow=TRUE,nrow=nrow(res$se),ncol=ncol(res$se))*sqrt(w)*dispersion
      }     
	}
	
  if(any(is.na(res$se)|(!is.finite(res$se)))){
    if(is.numeric(replaceNA.coeffWithin.se)){
      res$se[is.na(res$se)]=rep(as.double(replaceNA.coeffWithin.se),length.out=sum(is.na(res$se)))
    } else if(replaceNA.coeffWithin.se=="coeffMeans"){
      for(i in which(apply(res$se,2,function(x) any(is.na(x))) )){ 
        res$se[is.na(res$se[,i]),i]=sqrt(mean(res$se[,i]^2,na.rm=TRUE))
      }
      covAve=apply(res$covs,c(2,3),mean,na.rm=TRUE)
      for(i in 1:nrow(res$covs)){ 
        getNA=which(is.na(res$covs[i,,]))
        if(length(getNA)>0) {
          temp=res$covs[i,,]
          temp[getNA]=covAve[getNA]
          res$covs[i,,]=temp
        }
      }
    } else if(replaceNA.coeffWithin.se=="unitMeans") {
      for(i in which(apply(res$se,1,function(x) any(is.na(x))))) 
        res$se[i,is.na(res$se[i,])]=sqrt(mean(res$se[i,]^2,na.rm=TRUE))
    }
  }
  if(is.null(res$Z)) res$Z=matrix(,nrow(res$X),0)
	res
}

#############################
.testedWithin <- function(data,includeOnly=NULL, exclude=NULL){
  
  if(!is.null(includeOnly)){  
    data$coeffWithin = data$coeffWithin[,includeOnly]
    data$se = data$se[,includeOnly]
    data$df.mod = data$df.mod[,includeOnly]
    data$covs = data$covs[,includeOnly,includeOnly]
#     data$df.res <- data$df.res[,includeOnly,drop=FALSE]
#     data$dispersion <- data$dispersion[,includeOnly,drop=FALSE]
  }
  if(!is.null(exclude)){
    if(is.character(exclude))
      exclude=unlist(sapply(exclude,function(x) which(names(data$coeffWithin)==x)))
    if(length(exclude)>0){
      data$coeffWithin = data$coeffWithin[,-exclude]
      data$se = data$se[,-exclude]
      data$df.mod = data$df.mod[,-exclude]
      data$covs = data$covs[,-exclude,-exclude]
#       data$df.res <- data$df.res[,-exclude,drop=FALSE]
#       data$dispersion <- data$dispersion[,-exclude,drop=FALSE]
    }
  }  
  data
}

#############################
.getEmptyFlipMix <- function(colNames,idClust,ncoef,nUnit){
	if(missing(ncoef)) ncoef=length(colNames)
	if(missing(nUnit)) nUnit=length(unique(idClust))
    
	res=list(coeffWithin=matrix(NA,nUnit,ncoef,dimnames=list(as.character(idClust),colNames)))
	dispersion=matrix(NA,nUnit,min(ncoef,1),dimnames=list(as.character(idClust)))
	covs=array(NA,c(nUnit,ncoef,ncoef))
	dimnames(covs)[[1]]=idClust
	dimnames(covs)[[2]]<-dimnames(covs)[[3]]<-colNames
	res=list(coeffWithin=res[[1]],se=res$coeffWithin,df.mod=res$coeffWithin,df.res=dispersion,covs=covs,dispersion=dispersion)
	res
}

##############################
.unlist.obs2coeffWithin <- function(.out){
	#colNames = unlist(lapply(names(.out),function(modName) paste(modName,names(coefficients(.out[[modName]])),sep=".")))
	res=.getEmptyFlipMix(c(),idClust=rownames(.out[[1]]$coeffWithin))
	for(i in 1:length(.out)){
		colnames(.out[[i]]$coeffWithin)=paste(names(.out)[i],colnames(.out[[i]]$coeffWithin),sep=".")
		res$coeffWithin <- cbind(res$coeffWithin,.out[[i]]$coeffWithin)
		res$se <- cbind(res$se,.out[[i]]$se)
		res$df.mod <- cbind(res$df.mod,.out[[i]]$df.mod)
		}
	colnames(res$se) <- colnames(res$df.mod) <- colnames(res$coeffWithin)
	res$df.res=NULL
	res$covs=NULL
	res$dispersion=NULL
	res
}
