############################
# permutation Test of dependence 
# Y is the Nxp matrix of responces
# X is the Nxq matrix of predictors
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of permutations (IDs) be returned?
############################

.dependence.nptest <- function(data, perms=5000,  tail = NULL,statTest="t",separatedX=TRUE,testType="permutation",return.permIDs=TRUE,...){
	
	if(is.null(statTest)) 
		if(is.null(separatedX) || separatedX) 
			statTest="t" else
			statTest="F"

	if((testType=="rotation") && statTest%in%c("Fisher","Wilcoxon","Kruskal-Wallis","rank","chisq","chisq.separated") ){
	  warning("Rotations are not allowed for Fisher exact test, permutations will be used instead.")
	  testType="permutation"
	}

    if(statTest=="NA"){ #for missing values
	    #same function for permutation and rotation tests
	      test <- .NA.dependence.nptest
	  } else if(statTest%in%c("sum","t","coeff")){
				test <- .t.dependence.nptest		
		} else if(statTest=="F"){ #ANOVAtype test, 1 column for each column of Y summarizing the dependence with all Xs			
				test <- .F.dependence.nptest
		} else if(statTest=="Fisher"){
				test <- .fisher.dependence.nptest
		} else if(statTest%in%c("Wilcoxon","Kruskal-Wallis","rank")){
				test <- .rank.dependence.nptest
		} else if(statTest%in%c("chisq","chisq.separated")){
			test <- .chisq.dependence.nptest
# 		} else if(statTest%in%"Kolmogorov-Smirnov"){
# 			test <- .kolmogorov.dependence.nptest
		} else	{stop("This test statistic is not valid, nothing done."); return()}
  
if(statTest%in%c("Fisher","Wilcoxon","Kruskal-Wallis","rank","chisq","Kolmogorov-Smirnov")) {
			if(length(unique(data$Z))>1 ) warning("Covariates Z can not be used in this test. Use strata instread.")
	}
  environment(test) <- sys.frame(sys.nframe())
  out <- sys.frame(sys.nframe())
  return(out)
}

###########################################

.t.dependence.nptest <- function(){
  data <- .orthoZ(data) #if Z is.null, it doesn't make anything
	N=nrow(data$Y)
  perms <- make.permSpace(1:N,perms,return.permIDs=TRUE,testType=testType, Strata=data$Strata)
	#search for intercept in the model
  intercept=.getIntercept(data$X)
  if(statTest=="coeff") {
    data$X=data$X%*%solve(t(data$X)%*%data$X)
  }
  
	if(any(intercept)) {
		data$intercept=TRUE
		data$X=data$X[,!intercept,drop=FALSE]
	}
  
    
		#data$X=scale(data$X,scale=FALSE)
  permT=.prod.perms(data,perms,testType=testType)
	
	if(statTest%in%c("sum","coeff") )
		permT= .prod2sum(permT,data)
	else if(statTest=="t") {
		permT= .prod2t(permT,data)
		if(any(intercept) ) permT =permT * sqrt((N-1-sum(intercept))/(N-1))
		}

	colnames(permT) = .getTNames(data$Y,data$X,permT=permT)
	rownames(permT)=.getTRowNames(permT)		
	
	if(statTest=="sum") {
		center=as.vector(outer(colSums(data$X),colSums(data$Y),"*"))/N
		names(center)=colnames(permT)
		attributes(tail)$center=center
	}
	res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test=statTest))
}
			
.F.dependence.nptest <- function(){
	data <- .orthoZ(data)
	N=nrow(data$Y)	
	perms <- make.permSpace(1:N,perms,testType=testType,return.permIDs=TRUE, Strata=data$Strata)
	#search for intercept in the model
	intercept=.getIntercept(data$X) 
	if(any(intercept)) {
		data$X=data$X[,!intercept,drop=FALSE]
		data$X=scale(data$X,scale=FALSE)
		data$Y=scale(data$Y,scale=FALSE)
	}
	
	
	
#	environment(.prod2F) <- sys.frame(sys.nframe())
	P=.get.eigenv.proj.mat(data)
	permT=.prod.perms.P(data,perms,testType=testType,P)
	permT=.prod2F(permT,data)
  if(any(intercept)) {
    q=ncol(data$X)#-sum(intercept) #just a trick, it should be subtracetd to N, not to q
    permT =permT * (N-q-sum(intercept))/(N-q)
  }
	
	colnames(permT) = .getTNames(data$Y,permT=permT)
	rownames(permT)=.getTRowNames(permT)		
	res=list(permT=permT,perms=perms,tail=1,extraInfoPre=list(Test="F"))
}

#################################
.rank.dependence.nptest <- function(){
	data$Y=apply(data$Y,2,rank)
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE,testType=testType, Strata=data$Strata)
	
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) data$X=data$X[,!intercept,drop=FALSE]

	if(statTest=="rank")
			statTest=ifelse(ncol(data$X)>1,"Kruskal-Wallis","Wilcoxon")

	if(statTest=="Wilcoxon"){
		permT=.prod.perms(data,perms)
		permT=scale(permT,center=rep(colSums(data$X)*(nrow(data$X)+1)/2,ncol(data$Y)))
		colnames(permT)=.getTNames(data$Y,data$X)
		rownames(permT)=.getTRowNames(permT)
		res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="Wilcoxon"))
	} else { #kruskall-wallis (anova)
		data$X=scale(data$X,scale=FALSE)
		data$Y=scale(data$Y,scale=FALSE)
		P=.get.eigenv.proj.mat(data)
		permT=.prod.perms.P(data,perms,testType=testType,P)
		permT=.prod2F(permT,data)
		permT=permT *ncol(data$X)/(nrow(data$Y)-ncol(data$X)) 
		permT = (1+ permT^-1)^-1 *(nrow(data$Y)-1)
		res=list(permT=permT,perms=perms,tail=1,extraInfoPre=list(Test="Kruskal-Wallis")) 
		res
	}
	
}
#################################
.fisher.dependence.nptest <- function(){
	if(any(apply(data$Y,2,function(y) length(unique(y))!=2 ))) stop("Only factor or dichotomous variables are allowed with Fisher exact test. Nothing done.")
	
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE,testType=testType, Strata=data$Strata)
	
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) data$X=data$X[,!intercept,drop=FALSE]

#	data$X=scale(data$X,scale=FALSE)
	permT=.prod.perms(data,perms)
	
	colnames(permT)=.getTNames(data$Y)
	rownames(permT)=.getTRowNames(permT)		
	center=as.vector(colSums(data$X)%*%t(colSums(data$Y)))/N
	names(center)=colnames(permT)
	attributes(tail)<-list(center=center)
	res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="Fisher"))
}
################################
# .kolmogorov.nptest<-function(){
#   
#   }
####################################
.chisq.dependence.nptest<-function(){ 
	if(any(apply(data$Y,2,function(y) length(unique(y))!=2 ))) stop("Only factor or dichotomous variables are allowed with Chi Squared test. Nothing done.")
	
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE,testType=testType, Strata=data$Strata)
	
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) { 
		attrsX=attributes(data$X)[c("assign","factors")]
		data$X=data$X[,!intercept,drop=FALSE]
		attrsX$assign=attrsX$assign[attrsX$assign>0]
		attributes(data$X)[c("assign","factors")] <- attrsX
	}
#	data$X=scale(data$X,scale=FALSE)
	permT=.prod.perms(data,perms)
	expected=as.vector(colSums(data$X)%*%t(colSums(data$Y)))/N
	permT=scale(permT,center=expected,scale=sqrt(expected))
	if(statTest=="chisq.separated"){
		colnames(permT)=.getTNames(data$Y,data$X) 
	} else {
		newNames=apply(expand.grid(attributes(data$X)$factors,attributes(data$Y)$factors),1,paste,collapse="_|_")
		permT=sapply(unique(attributes(data$Y)$assign),
				function(idy) {
					whichCol= as.vector(outer(1:ncol(data$X),(which(attributes(data$Y)$assign==idy)-1)*ncol(data$X),"+"))
					rowSums(permT[,whichCol,drop=FALSE]^2)
					})
		colnames(permT)=newNames
	}
	rownames(permT)=.getTRowNames(permT)	
	if(statTest=="chisq")
	  res=list(permT=permT,perms=perms,tail= 1,extraInfoPre=list(Test="Chi Squared")) else
	res=list(permT=permT,perms=perms,tail= tail,extraInfoPre=list(Test="signed-Chisq"))
}

############# NA

.NA.dependence.nptest <- function(){
  #data <- .orthoZ(data) #if Z is.null, it doesn't make anything
  N=nrow(data$Y)
  perms <- make.permSpace(1:N,perms,return.permIDs=return.permIDs,testType=testType, Strata=data$Strata)
  b=NULL  
  data$X=data$X[,which(!.getIntercept(data$X)),drop=FALSE]
  
  data$Y=apply(data$Y,2,scale,scale=FALSE)
  
  if(ncol(data$X)>1) { #multiple predictors
    #browser()
      notNA=!is.na(data$Y)
      tObs= as.vector(unlist(sapply(1:ncol(data$Y),function(i){
        x=apply(data$X[notNA[,i],,drop=FALSE],2,scale,scale=FALSE)
        ei=eigen(t(x)%*%x)
        if(any(ei$values<=1E-10))
          NA else { #at least a constant term in X
          F=sum(((t(data$Y[notNA[,i],i])%*%x%*% ei$vectors) ^2) %*% (ei$values^-1))
          #equivalent to t(data$Y[notNA[,i],i])%*%x%*% (ei$vectors %*% diag(ei$values^-1) %*% t(ei$vectors) ) %*%t(x)%*%data$Y[notNA[,i],i] #explained dev
          #but faster
        } 
      })))
      rm(notNA)

      permT=matrix(,perms$B,length(tObs))
      permT[1,]=tObs
      for(b  in (1:(perms$B-1))) {
                    Yperm=perms$rotFunct(b)
                    notNA=!is.na(Yperm)
                    permT[b+1,]=as.vector(unlist(sapply(1:ncol(data$Y),function(i) {
                      if(is.na(tObs[i])) NA else 
                        {x=apply(data$X[notNA[,i],,drop=FALSE],2,scale,scale=FALSE)
                         ei=eigen(t(x)%*%x)
                         if(any(ei$values<=1E-10)) {#browser();
                           NA} else { #at least a constant term in X
                             F=sum(((t(Yperm[notNA[,i],i])%*%x%*% ei$vectors) ^2) %*% (ei$values^-1))
                             # equivalent to : t(Yperm[notNA[,i],i])%*%x%*% (ei$vectors %*% diag(ei$values^-1) %*% t(ei$vectors) ) %*%t(x)%*%Yperm[notNA[,i],i] #explained dev
                           }
                      }
                    })))
                  }
      rm(tObs)
      permT=t( t(permT)/(apply(data$Y,2,var,na.rm=TRUE)*(apply(notNA,2,sum)-1) -t(permT)) / ( (ncol(data$X)) / (apply(notNA,2,sum) - ncol(data$X)-1) )  )
      
      colnames(permT) = .getTNames(data$Y,permT=permT)
  } else { #only one predictor
    data$Y=scale(data$Y, center = TRUE, scale = FALSE)
      notNA=!is.na(data$Y)
      tObs= as.vector(unlist(sapply(1:ncol(data$Y),function(i) {
        if(length(unique(data$X[notNA[,i],]))==1) NA  #X is constant
        else t(scale(data$X[notNA[,i],]))%*%data$Y[notNA[,i],i]
      })))
      rm(notNA)
      permT=matrix(,perms$B,length(tObs))
      permT[1,]=tObs
      for(b in (1:(perms$B-1))){
                Yperm=perms$rotFunct(b)
                notNA=!is.na(Yperm)
                permT[b+1,]=as.vector(unlist(sapply(1:ncol(data$Y),function(i){ 
                  if(is.na(tObs[i]) || (length(unique(data$X[notNA[,i],]))==1) ) 
                    NA #X is constant
                  else t(scale(data$X[notNA[,i],]))%*%Yperm[notNA[,i],i] 
                })))
              }
      rm(tObs)
      permT <- t(t(permT)*sqrt(apply(notNA,2,sum)-1))
      #browser() #want to transform to t stat
      #permT <- t( t(permT)/sqrt((apply(data$Y,2,var,na.rm=TRUE)*(apply(notNA,2,sum)-1) -t(permT)) / 
      #                            (apply(notNA,2,sum) - 2) ) )  
      
        
  colnames(permT) = .getTNames(data$Y,data$X,permT=permT)
  }
  rownames(permT)=.getTRowNames(permT)  	
  res=list(permT=permT,perms=perms,tail=tail,
           extraInfoPre=list(Test=ifelse(ncol(data$X)>1,"F-NA","uni-NA"),NValidPerms=apply(!is.na(permT),2,sum)))
}
