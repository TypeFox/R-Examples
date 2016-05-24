"predict.MCMCglmm"<-function(object, newdata=NULL, marginal=object$Random$formula, type="response", interval="none", level=0.95, it=NULL, posterior="all", verbose=FALSE, ...){

 # warning("predict.MCMCglmm is still developmental - be careful")

  rm.obs<-c()

  if(interval=="prediction"){
      post.pred<-t(simulate(object=object, nsim=nrow(object$Sol), newdata=newdata, marginal=marginal, type=type, it=it, posterior=posterior, verbose=verbose))

  }else{

    if(type%in%c("response", "terms")==FALSE){stop("type must be response or terms")}
    if(interval%in%c("none", "confidence", "prediction")==FALSE){stop("interval must be none, confidence or prediction")}

    if(!is.null(posterior)){
      if(!posterior%in%c("distribution", "mean", "mode", "all")){
        stop("posterior argument must be either distribution, mean, mode or all")
      }
    }

    if(!is.null(marginal)){
      if(class(marginal)!="formula"){stop("marginal should be NULL or a formula")}
    }

    if(!is.null(it)){
      if(length(it)>1){stop("it should be an integer")}
      if(it>nrow(object$X) | it<1){stop("it should be less than or equal to the number of iterations")}
    }

    rcomponents<-split.direct.sum(as.character(object$Random$formula)[2])
    mcomponents<-split.direct.sum(as.character(marginal)[2])

    if(length(rcomponents)!=length(object$Random$nrt)){stop("sorry - not implented for covu models")}
    if(any(mcomponents%in%rcomponents==FALSE)){stop("marginal formula does not correspond to model formula")}

    marginalise<-rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt)

    if(!is.null(newdata)){
      suppressWarnings(object2<-MCMCglmm(fixed=object$Fixed$formula, random=object$Random$formula, rcov=object$Residual$formula, family=object$Residual$original.family, data=newdata, nitt=1, thin=1, burnin=0, ginverse=object$ginverse, verbose=FALSE, pr=any(marginalise==0)))
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
      object2$Sol<-object$Sol[,c(find.fixed, find.random)]
      find.vcv<-match(colnames(object2$VCV), colnames(object$VCV))
      if(any(is.na(find.vcv))){stop("model for newdata has (co)variance terms not in original model")}
      object2$VCV<-object$VCV[,find.vcv,drop=FALSE]
      if(!is.null(object2$CP)){
        find.cp<-match(colnames(object2$CP), colnames(object$CP))
        if(any(is.na(find.cp))){stop("model for newdata has cutpoints not in original model")}
        object2$CP<-object$CP[,find.cp,drop=FALSE]
      }
      object<-object2
      rm(object2)
    }

    if(posterior=="mean"){
      object$VCV<-matrix(colMeans(object$VCV), 1, ncol(object$VCV))
      it<-1
    }
    if(posterior=="mode"){
      object$VCV<-matrix(posterior.mode(object$VCV, ...), 1, ncol(object$VCV))
      it<-1
    }
    if(is.null(it)){
      if(posterior=="distribution"){
        it<-sample(nrow(object$Sol), 1)
      }
      if(posterior=="all"){
        it<-1:nrow(object$Sol)
      }
    }
    object$Sol<-object$Sol[it,,drop=FALSE]
    object$VCV<-object$VCV[it,,drop=FALSE]
    if(!is.null(object$Lambda)){
      object$Lambda<-object$Lambda[it,,drop=FALSE]
    }
    if(!is.null(object$CP)){
      object$CP<-object$CP[it,,drop=FALSE]
    }

    if(is.null(object$X)){
      stop("fixed effect design matrix not saved: pass saveX=TRUE to MCMCglmm")
    }

    if(any(marginalise==0) & dim(object$Sol)[2]==dim(object$X)[2]){
      stop("posterior distribution of random effects not saved: pass pr=TRUE to MCMCglmm")
    }
    if(any(marginalise==0) & is.null(object$Z)){
      stop("random effect design matrix not saved: pass saveZ=TRUE to MCMCglmm")
    }


    if(is.null(object$Random$nfl)==FALSE){  # there are random effects
      st<-c(1,cumsum(object$Random$nrl*object$Random$nfl)+1)  # starting column for random effects of each component 
      st<-st[-length(st)]
      end<-cumsum(object$Random$nrl*object$Random$nfl)        # ennding column for random effects of each component 
      keep<-unlist(mapply(st[which(marginalise==0)], end[which(marginalise==0)], FUN=":"))    # random effects to be kept
    }else{
      keep<-NULL
    }

    if(!is.null(newdata) & !is.null(object$Random$nfl)){
      missing.random<-which(is.na(object$Sol[1,]))
      if(posterior=="mean" | posterior=="mode"){
        if(any(is.na(object$Sol))){
          object$Sol[missing.random]<-0
        }
      }else{
        dv<-1:ncol(object$Sol)
        cnt<-0
        for(j in 1:length(object$Random$nfl)){
          nfl<-object$Random$nfl[j]
          nrl<-object$Random$nrl[j]
          dv[object$Fixed$nfl+cnt+1:(nfl*nrl)]<-rep(diag(matrix(1:(nfl^2),nfl,nfl)), each=nrl)
          cnt<-cnt+(nfl*nrl)
        }
        if(any(is.na(object$Sol))){
          object$Sol[,missing.random,drop=FALSE]<-rnorm(nrow(object$Sol)*length(missing.random), 0, sqrt(object$VCV[,dv[missing.random]]))
        }
      }
    }

    object$Sol<-object$Sol[,c(1:object$Fixed$nfl, object$Fixed$nfl+keep),drop=FALSE]

    W<-cBind(object$X, object$Z)
    W<-W[,c(1:object$Fixed$nfl, object$Fixed$nfl+keep), drop=FALSE]
  
    post.pred<-t(apply(object$Sol,1, function(x){(W%*%x)@x}))

    if(type=="response"){
 
      if(any(object$family!="gaussian" & object$family!="cengaussian")){
         post.var<-buildV(object, marginal=marginal, diag=TRUE, it=NULL, posterior="all", ...)
      }

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

      normal.evd<-function(mu, v){
        int.foo<-function(x, mu, v){exp(-exp(x))*dnorm(x, mu, sqrt(v))}
        integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v)[[1]]
      }
      normal.zt<-function(mu, v){
        int.foo<-function(x, mu, v){exp(x)/(1-exp(-exp(x)))*dnorm(x, mu, sqrt(v))}
        integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v)[[1]]
      }  
      normal.logistic<-function(mu, v){
        int.foo<-function(x, mu, v){plogis(x)*dnorm(x, mu, sqrt(v))}
        integrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu,v)[[1]]
      } 
      normal.multilogistic<-function(mu, v){
        int.foo<-function(x, mu, v, i){(exp(x[i])/(1+sum(exp(x))))*prod(dnorm(x, mu, sqrt(v)))}
        res<-1:length(mu)
        for(i in 1:length(mu)){
           res[i]<-cubature::adaptIntegrate(int.foo, qnorm(0.0001, mu,sqrt(v)), qnorm(0.9999, mu,sqrt(v)), mu=mu, v=v, i=i)[[1]]
        }
        res
      } 

      if(!is.null(object$Lambda)){
        if(any(object$family!="gaussian")){stop("sorry - prediction not yet available for sir models applied to non-Gaussian data")}
        # bear in mind that if this is implemented the post.var is quadratic in the structural parameters
        I<-Diagonal(nrow(object$XL))
        post.pred<-t(sapply(1:nrow(object$Sol), function(i){as.vector(solve(I-object$XL%*%kronecker(object$Lambda[i,], I), post.pred[i,]))}))
      }

      if(any(object$family%in%c("poisson","cenpoisson"))){
        keep<-which(object$family%in%c("poisson","cenpoisson"))
        post.pred[,keep]<-exp(post.pred[,keep]+0.5*post.var[,keep])
      }

      if(any(object$family%in%c("ztpoisson"))){
        keep<-which(object$family%in%c("ztpoisson"))
        post.pred[,keep]<-mapply(post.pred[,keep], post.var[,keep], FUN=function(mu,v){normal.zt(mu,v)})
      }

      if(any(object$family%in%c("exponential","cenexponential","geometric","cengeometric"))){
        keep<-which(object$family%in%c("exponential","cenexponential","geometric","cengeometric"))
        post.pred[,keep]<-exp(-post.pred[,keep]+0.5*post.var[,keep])
      }

      if(any(object$family%in%c("ordinal"))){      

        nord<-unique(object$error.term[which(object$family=="ordinal")])
        cp.names<-substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", colnames(object$CP))-1)
        cp.names<-match(cp.names,unique(cp.names))

        for(k in 1:length(nord)){

          keep<-which(object$family%in%c("ordinal") & object$error.term==nord[k])
          CP<-cbind(-Inf, 0, object$CP[,which(cp.names==k), drop=FALSE], Inf)
          q<-matrix(0,dim(post.pred)[1], length(keep))

          for(i in 2:(dim(CP)[2]-1)){
            q<-q+(pnorm(CP[,i+1]-post.pred[,keep],0,sqrt(post.var[,keep]+1))-pnorm(CP[,i]-post.pred[,keep],0,sqrt(post.var[,keep]+1)))*(i-1)
          }

          post.pred[,keep]<-q
          rm(q)
        } 
      }

      if(any(object$family%in%c("threshold"))){      
      
        nord<-unique(object$error.term[which(object$family=="threshold")])
        cp.names<-substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", colnames(object$CP))-1)
        cp.names<-match(cp.names,unique(cp.names))

        for(k in 1:length(nord)){
          keep<-which(object$family%in%c("threshold") & object$error.term==nord[k])
          CP<-cbind(-Inf, 0, object$CP[,which(cp.names==k), drop=FALSE], Inf)
          q<-matrix(0,dim(post.pred)[1], length(keep))

          for(i in 2:(dim(CP)[2]-1)){
            q<-q+(pnorm(CP[,i+1]-post.pred[,keep],0,sqrt(post.var[,keep]))-pnorm(CP[,i]-post.pred[,keep],0,sqrt(post.var[,keep])))*(i-1)
          }

          post.pred[,keep]<-q
          rm(q)
        } 
      }

      for(k in unique(super.trait)){
        if(any(grepl("multinomial", object$Residual$family))){
          keep<-which(object$error.term%in%which(super.trait==k))
          size<-as.numeric(substr(object$family[keep], 12, nchar(object$family[keep])))
          for(j in 1:nrow(post.pred)){
            prob<-matrix(post.pred[j,keep], length(keep)/sum(super.trait==k), sum(super.trait==k))
            pvar<-matrix(post.var[j,keep], length(keep)/sum(super.trait==k), sum(super.trait==k))
            post.pred[j,keep]<-t(sapply(1:nrow(prob), function(x){normal.multilogistic(prob[x,], pvar[x,])}))*size
          }
        }

        if(any(grepl("hupoisson", object$Residual$family))){
          keep<-which(object$error.term%in%which(super.trait==k))
          keep<-keep[-c(1:(length(keep)/2))]
          rm.obs<-c(rm.obs, keep)
          post.pred[,keep]<-mapply(post.pred[,keep], post.var[,keep], FUN=function(mu,v){1-normal.logistic(mu,v)})
          post.pred[,keep-length(keep)]<-post.pred[,keep]*mapply(post.pred[,keep-length(keep)], post.var[,keep-length(keep)], FUN=function(mu,v){normal.zt(mu,v)})
        }
        if(any(grepl("zapoisson", object$Residual$family))){
          keep<-which(object$error.term%in%which(super.trait==k))
          keep<-keep[-c(1:(length(keep)/2))]
          rm.obs<-c(rm.obs, keep)
          post.pred[,keep]<-1-mapply(post.pred[,keep], post.var[,keep], FUN=function(mu,v){normal.evd(mu,v)})
          post.pred[,keep-length(keep)]<-post.pred[,keep]*mapply(post.pred[,keep-length(keep)], post.var[,keep-length(keep)], FUN=function(mu,v){normal.zt(mu,v)})
        }

        if(any(grepl("zipoisson", object$Residual$family[which(super.trait==k)]))){
          keep<-which(object$error.term%in%which(super.trait==k))
          keep<-keep[-c(1:(length(keep)/2))]
          rm.obs<-c(rm.obs, keep)
          post.pred[,keep]<-mapply(post.pred[,keep], post.var[,keep], FUN=function(mu,v){1-normal.logistic(mu,v)})
          post.pred[,keep-length(keep)]<-post.pred[,keep]*exp(post.pred[,keep-length(keep)]+post.var[,keep-length(keep)]/2)
        }
        if(any(grepl("zibinomial", object$Residual$family))){
          keep<-which(object$error.term%in%which(super.trait==k))
          keep<-keep[-c(1:(length(keep)/2))]
          rm.obs<-c(rm.obs, keep)
          post.pred[,keep]<-mapply(post.pred[,keep], post.var[,keep], FUN=function(mu,v){1-normal.logistic(mu,v)})
          post.pred[,keep-length(keep)]<-post.pred[,keep]*mapply(post.pred[,keep-length(keep)], post.var[,keep-length(keep)], FUN=function(mu,v){normal.logistic(mu,v)})
        }
      }
    }
  }
  
  if(length(rm.obs)>0){
    post.pred<-post.pred[,-rm.obs]
  }

  if(is.matrix(post.pred)){
    pred<-matrix(colMeans(post.pred), dim(post.pred)[2],1)
  }else{
    pred<-matrix(post.pred, length(post.pred),1)
  }  

  if(interval!="none"){
    pred<-cbind(pred, coda::HPDinterval(mcmc(post.pred), prob=level))   
    colnames(pred)<-c("fit", "lwr", "upr")
  }
  
  rownames(pred)<-1:dim(pred)[1]

  return(pred)
}  
