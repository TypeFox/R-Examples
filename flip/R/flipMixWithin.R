
.estimateSuSymmetry <- function (Y,S,n) {
  # browser()
  if(length(dim(S))>=3) S<-t(apply(S,1,diag))
  Su <- pmax(0,(colSums(Y^2)-colSums(S))/n)
}

#####################
.reff.symmetry.nptest <- function(data, perms=5000, statTest="Tnaive",  tail = NULL,...){
  test <- function(){
    n=colSums(data$Y!=0)
    p=ncol(data$Y)
    nmax=nrow(data$Y)
    perms <- make.signSpace(nmax,perms)
    #perms <- make.permSpace(nmax,perms,return.permIDs=FALSE,testType=rotation")
    if(("TH0est"%in%statTest) || ("TH1est"%in%statTest) || ("TBTWest"%in%statTest) ) {
      em2=colSums(data$Y^2)
      se2=data$se^2
      Sw=colSums(se2)
    }
    
    ### no estimates
    if("Tnaive"%in%statTest){
      permT <- rbind(rep(1,perms$n),perms$permID) %*% (data$Y%*%diag(1/n))
      permT = rbind(permT,-permT[nrow(permT):1,,drop=FALSE])
      rownames(permT)=.getTRowNames(permT)
      colnames(permT) = .getTNames(data$Y,permT=permT)
      colnames(permT)=paste(colnames(permT),"_Tnaive",sep="")
    } else permT=NULL
    
    ### H0 estimates
    if("TH0est"%in%statTest){
      Su <- .estimateSuSymmetry (data$Y,se2,n)
      w=sqrt(se2+matrix(Su,nrow(data$Y),p,byrow=TRUE))
      Y=data$Y/w
      
      permT0 <- rbind(rep(1,perms$n),perms$permID) %*% Y
      permT0=permT0/sqrt(n)
      colnames(permT0) = .getTNames(data$Y)
      permT0 = rbind(permT0,-permT0[nrow(permT0):1,,drop=FALSE])
      rownames(permT0)=.getTRowNames(permT0)
      colnames(permT0) = .getTNames(data$Y,permT=permT0)
      colnames(permT0)=paste(colnames(permT0),"_H0est",sep="")
    } else permT0=NULL;
    
    ### H1 estimates    
    if("TH1est"%in%statTest){
      flips=rbind(rep(1,perms$n),perms$permID)
      permT1=matrix(,nrow(flips)*2,ncol(data$Y))
      if(is.null(permT)){
        permT1[1:nrow(flips),]=rbind(rep(1,perms$n),perms$permID) %*% (data$Y%*%diag(1/n))
      } else permT1=permT
      
      for(ii in 1:nrow(flips)){
        W=sqrt(se2+pmax(0,matrix((em2-permT1[ii,]^2*n-Sw)/(n-1),nmax,p,byrow=TRUE)))%*%diag(sqrt(n))
        Y=data$Y/W
        permT1[ii,]=flips[ii,]%*%Y
      }
      colnames(permT1) = .getTNames(data$Y)
      permT1[nrow(permT1):(nrow(permT1)/2+1),] = -permT1[1:(nrow(permT1)/2),,drop=FALSE]
      rownames(permT1)=.getTRowNames(permT1)
      colnames(permT1) = .getTNames(data$Y,permT=permT1)
      colnames(permT1)=paste(colnames(permT1),"_H1est",sep="")
    } else permT1=NULL
    
    if("TBTWest"%in%statTest){
      flips=rbind(rep(1,perms$n),perms$permID)
      permTbtw=matrix(,nrow(flips)*2,ncol(data$Y))
      for(ii in nrow(flips):1){
        SuBtw=diag(.estimateSuMultiILS(Y=diag(flips[ii,])%*%data$Y,Z=data$X, S=data$covs))
        W=sqrt(se2+matrix(SuBtw,nmax,p,byrow=TRUE))%*%diag(sqrt(n))
        Y=data$Y/W
        permTbtw[ii,]=flips[ii,]%*%Y
      }
       colnames(permTbtw) = .getTNames(data$Y)
       permTbtw[nrow(permTbtw):(nrow(permTbtw)/2+1),] = -permTbtw[1:(nrow(permTbtw)/2),,drop=FALSE]
      rownames(permTbtw)=.getTRowNames(permTbtw)
      colnames(permTbtw) = .getTNames(data$Y,permT=permTbtw)
      colnames(permTbtw)=paste(colnames(permTbtw),"_BTWest",sep="")
    } else permTbtw=NULL
    
    return(list(permT=cbind(permT,permT0,permT1,permTbtw),
                perms=perms,tail=tail,
                extraInfoPre=list(est.Su=c(if("Tnaive"%in%statTest) rep(NA,ncol(permT)) else NULL,
                                        if("TH0est"%in%statTest) Su else NULL,
                                        if("TH1est"%in%statTest) pmax(0,(em2-colMeans(data$Y)^2*n-Sw)/(n-1)) else NULL,
                                        if("TBTWest"%in%statTest) SuBtw else NULL),
                                  Test=rep(statTest,each=p)
                                        )))
}
  
  environment(test) <- sys.frame(sys.nframe())  
  out <- sys.frame(sys.nframe())
}


######################################################

flipMixWithin <- function(modelWithin,units, X=~1, perms=1000, data=NULL, tail=NULL,
                          statTest=NULL,flipReturn,
                          Su=NULL, equal.se=FALSE,se=NA,replaceNA.coeffWithin=0,
                          replaceNA.coeffWithin.se=Inf, ...) {
  
  otherParams= list(...)
  statTest=match.arg(statTest,c("TH0est","Tnaive","TH1est","TBTWest"),several.ok=TRUE)
  if(missing(flipReturn)||is.null(flipReturn)) 
    flipReturn=list(permT=TRUE,permP=FALSE,permSpace=FALSE,data=FALSE)
  
  if(is.null(statTest) ) if(is.null(otherParams$separatedX)   || otherParams$separatedX)   { statTest="t" } else statTest="F"
  
  
  # store the call
call <- match.call()  
  if(!(is.list(data) && (!is.data.frame(data)))) {    
    data<-obs2coeffWithin(modelWithin,Z=NULL,X=X,units=units, data=data,equal.se=equal.se,se=se,
                          replaceNA.coeffWithin=replaceNA.coeffWithin,replaceNA.coeffWithin.se=replaceNA.coeffWithin.se,...)
    if(is.null(data$se)) data$se=-t(apply(data$covs,1,diag))
  }
  rm(modelWithin)
  #########
  N = nrow(data$coeffWithin)
  p = ncol(data$coeffWithin)
  
  names(data)[names(data)=="coeffWithin"]="Y"
  
  res=.reff.symmetry.nptest(data, perms=perms, statTest=statTest,  tail = tail,...)
  out=res$test()
  res=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=res)
}