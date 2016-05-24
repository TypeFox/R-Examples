"simulate.MCMCglmm"<-function(object, nsim = 1, seed = NULL, newdata=NULL, marginal=object$Random$formula, type="response", it=NULL, posterior="all", verbose=FALSE, ...){

  if(!is.null(seed)){
    if(!is.integer(seed)){
      stop("seed must be an integer")
    }else{
      set.seed(seed)
    }
  }

  if(is.null(object$X)){
    stop("fixed effect design matrix not saved: pass saveX=TRUE to MCMCglmm")
  }

  if(!is.null(it)){
    if(length(it)>1){stop("'it' should be an integer")}
    if(it>nrow(object$Sol) | it<1){stop("'it' should be less than or equal to the number of iterations")}
  }

  if(!is.null(marginal)){
    if(class(marginal)!="formula"){stop("marginal should be NULL or a formula")}
  }

  if(!is.null(posterior)){
    if(!posterior%in%c("distribution", "mean", "mode", "all")){
      stop("posterior argument must be either distribution, mean, mode or all")
    }
    if(posterior=="all" & nsim>nrow(object$Sol)){
       stop("nsim must be less than the number of saved iterations if posterior='all'")
    }
  }

  if(type%in%c("response", "terms")==FALSE){stop("type must be response or terms")}

  rcomponents<-split.direct.sum(as.character(object$Random$formula)[2])
  mcomponents<-split.direct.sum(as.character(marginal)[2])

  if(length(rcomponents)!=length(object$Random$nrt)){stop("sorry - not implented for covu models")}
  if(any(mcomponents%in%rcomponents==FALSE)){stop("marginal formula does not correspond to model formula")}

  marginalise<-rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt)

  if(!is.null(newdata)){
   suppressWarnings(object2<-MCMCglmm(fixed=object$Fixed$formula, random=object$Random$formula, rcov=object$Residual$formula, family=object$Residual$original.family, data=newdata, nitt=1, thin=1, burnin=0, ginverse=object$ginverse, verbose=FALSE, pr=any(marginalise==0), start=list(QUASI=FALSE)))
   find.fixed<-match(colnames(object2$Sol)[1:object2$Fixed$nfl], colnames(object$Sol))
   find.random<-match(colnames(object2$Sol)[-c(1:object2$Fixed$nfl)], colnames(object$Sol))
   if(any(is.na(find.fixed))){stop("model for newdata has fixed effects not present in original model")}
   if(verbose){
     if(any(is.na(find.random))){
       missing.random<-colnames(object2$Sol)[which(is.na(find.random))+object2$Fixed$nfl]
       warning(paste("model for newdata has random effects not present in original model:", paste(missing.random, collapse=", ")))
     }
     missing.fixed<-which(!colnames(object$Sol)[1:object$Fixed$nfl]%in%colnames(object2$Sol)[1:object2$Fixed$nfl])
     if(length(missing.fixed)>0){
        missing.fixed<-colnames(object$Sol)[1:object$Fixed$nfl][missing.fixed]
        warning(paste("original model has fixed effects not present in newdata:", paste(missing.fixed, collapse=", ")))
     }
   }
   object2$Sol<-object$Sol[,c(find.fixed, find.random), drop=FALSE]
   find.vcv<-match(colnames(object2$VCV), colnames(object$VCV))
   if(any(is.na(find.vcv))){stop("model for newdata has (co)variance terms not in original model")}
   object2$VCV<-object$VCV[,find.vcv, drop=FALSE]
   if(!is.null(object2$CP)){
     find.cp<-match(colnames(object2$CP), colnames(object$CP))
     if(any(is.na(find.cp))){stop("model for newdata has cutpoints not in original model")}
     object2$CP<-object$CP
   }
   object<-object2
   rm(object2)
  }

  if(any(marginalise==0) & dim(object$Sol)[2]==dim(object$X)[2]){
     stop("posterior distribution of random effects not saved: pass pr=TRUE to MCMCglmm")
  }
  if(any(marginalise==0) & is.null(object$Z)){
     stop("random effect design matrix not saved: pass saveZ=TRUE to MCMCglmm")
  }

  if(!is.null(object$Random$nfl)){                                                              # there are random effects
    st<-c(1,cumsum(object$Random$nrl*object$Random$nfl)+1)                                      # starting column for random effects of each component 
    st<-st[-length(st)]
    end<-cumsum(object$Random$nrl*object$Random$nfl)                                            # ending column for random effects of each component 
    keep<-c(unlist(mapply(st[which(marginalise==0)], end[which(marginalise==0)], FUN=":")))     # random effects to be kept
  }else{
    keep<-NULL
  }

  if(!is.null(newdata) & !is.null(object$Random$nfl)){
    missing.random<-which(is.na(object$Sol[1,]))
    if(length(missing.random)>0){
      keep<-setdiff(keep,missing.random-object$Fixed$nfl)
    }
  }

  if(!is.null(it)){
      object$Sol<-object$Sol[it,,drop=FALSE]
      object$VCV<-object$VCV[it,,drop=FALSE]
      if(!is.null(object$Lambda)){
        object$Lambda<-object$Lambda[it,,drop=FALSE]
      }
  }else{
    if(posterior=="mean"){
      object$Sol<-matrix(colMeans(object$Sol), 1, ncol(object$Sol))
      if(any(is.na(object$Sol))){
         object$Sol[which(is.na(object$Sol))]<-0
      }
      object$VCV<-matrix(colMeans(object$VCV), 1, ncol(object$VCV))
      if(!is.null(object$Lambda)){
        object$Lambda<-matrix(colMeans(object$Lambda), 1, ncol(object$Lambda))
      }
    }

    if(posterior=="mode"){
      object$Sol<-matrix(posterior.mode(object$Sol, ...), 1, ncol(object$Sol))
      if(any(is.na(object$Sol))){
         object$Sol[which(is.na(object$Sol))]<-0
      }
      object$VCV<-matrix(posterior.mode(object$VCV, ...), 1, ncol(object$VCV))
      if(!is.null(object$Lambda)){
        object$Lambda<-matrix(posterior.mode(object$Lambda), 1, ncol(object$Lambda))
      }
    }
  }


  ynew<-matrix(NA, nrow(object$X), nsim)
  unew<-matrix(NA,sum(object$Random$nfl*object$Random$nrl),1)
  enew<-matrix(NA,sum(object$Residual$nfl*object$Residual$nrl),1)

  it<-sample(1:nrow(object$Sol), nsim, replace=posterior!="all")

  super.trait<-1:length(object$Residual$family)
  cnt<-1
  i<-1
  while(i<=length(super.trait)){
    if(!grepl("hu|zi|za|multinomial", object$Residual$family[i])){
      super.trait[i]<-cnt
      i<-i+1
      cnt<-cnt+1
    }else{
      if(grepl("multinomial", object$Residual$family[i])){
        nm<-as.numeric(substr(object$Residual$family[i], 12,nchar(object$Residual$family[i])))-1
        super.trait[i+1:nm-1]<-cnt
        i<-i+nm
        cnt<-cnt+1
      }else{
        super.trait[i]<-cnt
        super.trait[i+1]<-cnt
        i<-i+2
        cnt<-cnt+1
      }
    }
  }
    
  rm.obs<-c()

  for(i in 1:nsim){

    cnt<-0
    cnt2<-0

    if(!is.null(object$Random$nfl)){
       if(length(keep)!=ncol(object$Z)){
         for(j in 1:length(object$Random$nfl)){

           nfl<-object$Random$nfl[j]
           nrl<-object$Random$nrl[j]
           nat<-object$Random$nat[j]

           Y<-matrix(rnorm(nrl*nfl),nrl,nfl)

           if(nat==0){
             unew[1:(nrl*nfl)+cnt]<-as.vector(Y%*%chol(matrix(object$VCV[it[i],cnt2+1:(nfl^2)],nfl,nfl)))
           }else{
             unew[1:(nrl*nfl)+cnt]<-as.vector(solve(chol(forceSymmetric(object$ginverse[[nat]])), Y%*%chol(matrix(object$VCV[it[i],cnt2+1:(nfl^2)],nfl,nfl))))
           }
           cnt<-cnt+(nrl*nfl)
           cnt2<-cnt2+nfl^2
        }
      }else{
        cnt2<-sum(object$Random$nfl[1:length(object$Random$nfl)]^2)
      }
    }

    cnt<-0

    for(j in 1:length(object$Residual$nfl)){

      nfl<-object$Residual$nfl[j]
      nrl<-object$Residual$nrl[j]

      Y<-matrix(rnorm(nrl*nfl),nrl,nfl)

      enew[1:(nrl*nfl)+cnt]<-as.vector(Y%*%chol(matrix(object$VCV[it[i],cnt2+1:(nfl^2)],nfl,nfl)))

      cnt<-cnt+(nrl*nfl)
      cnt2<-cnt2+nfl^2
    }

    ynew[,i]<-as.vector(object$X%*%object$Sol[it[i],1:ncol(object$X)]+object$ZR%*%enew)

    if(!is.null(keep)){
        ynew[,i]<-ynew[,i]+as.vector(object$Z[,keep, drop=FALSE]%*%t(object$Sol[it[i],object$Fixed$nfl+keep, drop=FALSE]))
    }

    if(!is.null(object$Random$nfl)){
      if(length(keep)!=ncol(object$Z)){
        if(length(keep)==0){
          ynew[,i]<-ynew[,i]+as.vector(object$Z%*%unew)
        }else{

          ynew[,i]<-ynew[,i]+as.vector(object$Z[,-keep, drop=FALSE]%*%unew[-keep, drop=FALSE])
        }
      }
    }

    if(!is.null(object$Lambda)){
      ynew[,i]<-as.vector(solve(Diagonal(nrow(object$XL))-object$XL%*%kronecker(object$Lambda[it[i],], Diagonal(nrow(object$XL))), ynew[,i]))
    }

    if(type=="response"){

      if(any(object$family%in%c("poisson", "cenpoisson"))){
        trans<-which(object$family=="poisson" | object$family=="cenpoisson")
        ynew[trans,i]<-rpois(length(trans), exp(ynew[trans,i]))
      }
      if(any(object$family%in%c("exponential", "cenexponential"))){
        trans<-which(object$family=="exponential" | object$family=="cenexponential")
        ynew[trans,i]<-rexp(length(trans), exp(ynew[trans,i]))
      }
      if(any(object$family%in%c("geometric"))){
        trans<-which(object$family=="geometric")
        ynew[trans,i]<-rgeom(length(trans), plogis(ynew[trans,i]))
      }
      if(any(object$family%in%c("ztpoisson"))){
        trans<-which(object$family=="ztpoisson")
        ynew[trans,i]<-qpois(runif(length(trans), dpois(0, exp(ynew[trans,i])), 1), exp(ynew[trans,i]))  
      }
      if(any(object$family%in%"ordinal")){

        nord<-unique(object$error.term[which(object$family=="ordinal")])
        cp.names<-substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", colnames(object$CP))-1)
        cp.names<-match(cp.names,unique(cp.names))

        for(k in 1:length(nord)){
          trans<-which(object$family=="ordinal" & object$error.term==nord[k])
          CP<-c(-Inf, 0, object$CP[it[i],which(cp.names==k)], Inf)
          q<-matrix(NA,length(trans), length(CP)-1)
          for(j in 1:(length(CP)-1)){
            q[,j]<-pnorm(CP[j+1]-ynew[trans,i], 0, 1)-pnorm(CP[j]-ynew[trans,i], 0, 1)
          }
          ynew[trans,i]<-apply(q, 1, function(x){sample(1:ncol(q), 1, prob=x)})-1
        }
      }
      if(any(object$family%in%"threshold")){

        nord<-unique(object$error.term[which(object$family=="threshold")])
        cp.names<-substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", colnames(object$CP))-1)
        cp.names<-match(cp.names,unique(cp.names))

        for(k in 1:length(nord)){
          trans<-which(object$family=="threshold" & object$error.term==nord[k])
          ynew[trans,i]<-as.numeric(cut(ynew[trans,i], c(-Inf, 0, object$CP[it[i],which(cp.names==k), drop=FALSE], Inf)))-1         
        }
      }

      for(k in unique(super.trait)){
        if(any(grepl("multinomial", object$Residual$family[which(super.trait==k)]))){
          trans<-which(object$error.term%in%which(super.trait==k))
          prob<-matrix(ynew[trans,i], length(trans)/sum(super.trait==k), sum(super.trait==k))
          size<-as.numeric(substr(object$family[trans], 12, nchar(object$family[trans])))
          ynew[trans,i]<-t(sapply(1:nrow(prob), function(x){rmultinom(1, size=size[x], prob= c(1,exp(prob[x,]))/(1+sum(exp(prob[x,]))))}))[,-1]
        }
        if(any(grepl("hupoisson", object$Residual$family[which(super.trait==k)]))){
          trans<-which(object$error.term%in%which(super.trait==k))
          prob<-matrix(ynew[trans,i], length(trans)/sum(super.trait==k), sum(super.trait==k))
          if(i==1){
            rm.obs<-c(rm.obs, trans[-c(1:(length(trans)/sum(super.trait==k)))])
          }
          prob[,2]<-rbinom(nrow(prob), 1, 1-plogis(prob[,2]))
          prob[,2][which(prob[,2]==1)]<-qpois(runif(sum(prob[,2]==1), dpois(0, exp(prob[,1][which(prob[,2]==1)])), 1), exp(prob[,1][which(prob[,2]==1)]))
          ynew[trans,i]<-prob[,2]
        }
        if(any(grepl("zapoisson", object$Residual$family[which(super.trait==k)]))){
          trans<-which(object$error.term%in%which(super.trait==k))
          prob<-matrix(ynew[trans,i], length(trans)/sum(super.trait==k), sum(super.trait==k))
          if(i==1){
            rm.obs<-c(rm.obs, trans[-c(1:(length(trans)/sum(super.trait==k)))])
          }
          prob[,2]<-rbinom(nrow(prob), 1, pexp(exp(prob[,2])))
          prob[,2][which(prob[,2]==1)]<-qpois(runif(sum(prob[,2]==1), dpois(0, exp(prob[,1][which(prob[,2]==1)])), 1), exp(prob[,1][which(prob[,2]==1)]))
          ynew[trans,i]<-prob[,2]
        }
        if(any(grepl("zipoisson", object$Residual$family[which(super.trait==k)]))){
          trans<-which(object$error.term%in%which(super.trait==k))
          if(i==1){
            rm.obs<-c(rm.obs, trans[-c(1:(length(trans)/sum(super.trait==k)))])
          }
          prob<-matrix(ynew[trans,i], length(trans)/sum(super.trait==k), sum(super.trait==k))
          prob[,2]<-rbinom(nrow(prob), 1, 1-plogis(prob[,2]))
          prob[,2][which(prob[,2]==1)]<-rpois(sum(prob[,2]==1), exp(prob[,1][which(prob[,2]==1)]))
          ynew[trans,i]<-prob[,2]
        }
        if(any(grepl("zibinomial", object$Residual$family[which(super.trait==k)]))){
          trans<-which(object$error.term%in%which(super.trait==k))
          if(i==1){
            rm.obs<-c(rm.obs, trans[-c(1:(length(trans)/sum(super.trait==k)))])
          }
          prob<-matrix(ynew[trans,i], length(trans)/sum(super.trait==k), sum(super.trait==k))
          prob[,2]<-rbinom(nrow(prob), 1, 1-plogis(prob[,2]))
          prob[,2][which(prob[,2]==1)]<-rbinom(sum(prob[,2]==1), 1, plogis(prob[,1][which(prob[,2]==1)]))
          ynew[trans,i]<-prob[,2]
        }
      }
    }
  }
  if(length(rm.obs)>0){
    ynew<-ynew[-rm.obs,]
  }
  if(nsim==1){
    ynew<-as.vector(ynew)
  }
  return(ynew)
}  
