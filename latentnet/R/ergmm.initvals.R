ergmm.initvals <- function(model,user.start,prior,control){
  if(control[["verbose"]]) cat("Generating initial values for MCMC:\n")
  need.to.fit<-list(beta=model[["p"]]>0 && is.null(user.start[["beta"]]), ## beta
                    Z=model[["d"]]>0 && is.null(user.start[["Z"]]), ## Z
                    sender=model[["sender"]] && is.null(user.start[["sender"]]), ## sender
                    receiver=model[["receiver"]] && is.null(user.start[["receiver"]]), ## receiver
                    sociality=model[["sociality"]] && is.null(user.start[["sociality"]]),  ## sociality
                    Z.var=model[["d"]]>0 && is.null(user.start[["Z.var"]]),
                    Z.mean=model[["G"]]>0 && is.null(user.start[["Z.mean"]]),
                    Z.K=model[["G"]]>0 && is.null(user.start[["Z.K"]]),
                    Z.pK=model[["G"]]>0 && is.null(user.start[["pZ.K"]]),
                    sender.var=model[["sender"]] && is.null(user.start[["sender.var"]]), ## sender
                    receiver.var=model[["receiver"]] && is.null(user.start[["receiver.var"]]), ## receiver
                    sociality.var=model[["sociality"]] && is.null(user.start[["sociality.var"]]),
                    dispersion=model[["dispersion"]] && is.null(user.start[["dispersion"]])
                    )

  Yg<- model[["Yg"]]
  n <- network.size(Yg)
  p <- model[["p"]]
  d <- model[["d"]]
  G <- model[["G"]]

  Ym<-getYm(Yg,model[["response"]])
  
  Ym01<-Ym>mean(Ym,na.rm=TRUE)
  mode(Ym01)<-"numeric"
  
  pm<-user.start
  
  if(need.to.fit[["Z"]]){
    if(control[["verbose"]]) cat("Computing geodesic distances... ")
    D <- ergmm.geodesicmatrix(model)
    D[is.infinite(D)]<-2*n
    if(control[["verbose"]]) cat("Finished.\n")
    if(control[["verbose"]]) cat("Computing MDS locations... ")
    pm[["Z"]] <- cmdscale(D,model[["d"]])
    if(control[["verbose"]]) cat("Finished.\n")
  }

  if(control[["verbose"]]) cat("Computing other initial values... ")
  
  if("Z" %in% names(pm)) {
    i.keep<-mahalanobis(pm[["Z"]],0,cov(pm[["Z"]]))<20
    if(need.to.fit[["Z"]]) pm[["Z"]][!i.keep,]<-0
  }

  if(need.to.fit[["Z.K"]]){
    mbc1<-mbc.VII.EM(G,pm[["Z"]][i.keep,])
    pm[["Z.K"]]<-integer(n)
    pm[["Z.K"]][i.keep]<-mbc1[["Z.K"]]
    pm[["Z.K"]][!i.keep]<-which.max(tabulate(mbc1[["Z.K"]]))
  }
  
  if(need.to.fit[["Z.pK"]]){
    pm[["Z.pK"]]<-tabulate(pm[["Z.K"]])/n
  }
  
  if(need.to.fit[["Z.var"]]){
    if(!is.null(pm[["Z.K"]])) pm[["Z.var"]]<-sapply(1:G,function(g) var(c(subset(pm[["Z"]][i.keep,],pm[["Z.K"]][i.keep]==g))))
    else pm[["Z.var"]]<-var(c(pm[["Z"]][i.keep,]))
  }

  if(need.to.fit[["Z.mean"]]){
    pm[["Z.mean"]]<-do.call(rbind,lapply(1:G,function(g) apply(subset(pm[["Z"]][i.keep,,drop=FALSE],pm[["Z.K"]][i.keep]==g),2,mean)))
  }

  logit<-function(p) log(p/(1-p))
  
  if(need.to.fit[["beta"]]){
    if(model[["intercept"]])
      pm[["beta"]]<-logit(mean(Ym01,na.rm=TRUE))+if(!is.null(pm[["Z"]]))mean(as.matrix(dist(pm[["Z"]]))) else 0
    pm[["beta"]]<-c(pm[["beta"]],if(model[["intercept"]]) prior[["beta.mean"]][-1] else prior[["beta.mean"]])
  }

  bayes.prop<-function(x) (sum(x,na.rm=TRUE)+1)/(length(na.omit(x))+2)
  
  if(need.to.fit[["sociality"]]){
    pm[["sociality"]]<-logit((apply(Ym01,1,bayes.prop)+apply(Ym01,2,bayes.prop))/2)
    pm[["sociality"]]<-pm[["sociality"]]-mean(pm[["sociality"]])
  } else{
    if(need.to.fit[["sender"]]){
      pm[["sender"]]<-logit(apply(Ym01,1,bayes.prop))
      pm[["sender"]]<-pm[["sender"]]-mean(pm[["sender"]])
    }
    if(need.to.fit[["receiver"]]){
      pm[["receiver"]]<-logit(apply(Ym01,2,bayes.prop))
      pm[["receiver"]]<-pm[["receiver"]]-mean(pm[["receiver"]])
    }
  }

  if(need.to.fit[["sociality.var"]]){
    pm[["sociality.var"]]<-var(pm[["sociality"]])
  }

  if(need.to.fit[["sender.var"]]){
    pm[["sender.var"]]<-var(pm[["sender"]])
  }

  if(need.to.fit[["receiver.var"]]){
    pm[["receiver.var"]]<-var(pm[["receiver"]])
  }

  if(need.to.fit[["dispersion"]]){
    pm[["dispersion"]]<-1
  }

  if(control[["verbose"]]) cat("Finished.\n")
  
  if(control[["verbose"]]) cat("Finding the conditional posterior mode... ")
  if(control[["refine.user.start"]]){
    need.to.fit<-list(beta=model[["p"]]>0, ## beta
                      Z=model[["d"]]>0, ## Z
                      sender=model[["sender"]], ## sender
                      receiver=model[["receiver"]], ## receiver
                      sociality=model[["sociality"]],  ## sociality
                      Z.var=model[["d"]]>0,
                      Z.mean=model[["G"]]>0,
                      Z.K=model[["G"]]>0,
                      Z.pK=model[["G"]]>0,
                      sender.var=model[["sender"]], ## sender
                      receiver.var=model[["receiver"]], ## receiver
                      sociality.var=model[["sociality"]],
                      dispersion=model[["dispersion"]]
                      )
    user.start<-list()
  }
  for(i in 1:control[["mle.maxit"]]){
    if(control[["verbose"]]>1) cat(i,"")
    pm.old<-pm
    pm<-find.mpe(model,pm,
                 given=.merge.lists(list(Z.K=pm[["Z.K"]]),user.start),
                 prior=prior,control=control,fit.vars=need.to.fit)
    if(is.null(pm)) stop("Problem fitting. Starting values may have to be supplied by the user.")
    if(need.to.fit[["Z.K"]])pm[["Z.K"]]<-try(mbc.VII.EM(G,pm[["Z"]],resume=list(Z.mean=pm[["Z.mean"]],Z.var=pm[["Z.var"]],Z.pK=pm[["Z.pK"]]))[["Z.K"]])
    if(inherits(pm[["Z.K"]],"try-error")) stop("Unable to find an initial clustering. Try fitting a model with fewer clusters, or specifying initial clusters manually.")
    if(isTRUE(all.equal(pm.old,pm))) break
  }
  if(control[["verbose"]]) cat("Finished.\n")

  pm
}

mbc.VII.EM<-function(G,Z,EM.maxit=200,EM.tol=.Machine$double.eps^0.5,EM.maxstarts=15,resume=NULL){
  Z<-cbind(Z)
  n<-dim(Z)[1]
  d<-dim(Z)[2]

  for(attempt in 1:EM.maxstarts){
    if(is.null(resume) || attempt>2){
      if(G>1){
        cl<-kmeans(Z,G,nstart=EM.maxstarts)
        theta<-list(Z.mean = cl[["centers"]],
                    Z.var = cl[["withinss"]]/cl[["size"]],
                    Z.pK = cl[["size"]]/sum(cl[["size"]]))
        ## If the initial hard clustering produces a singleton
        ## cluster, boost the variance to avoid degeneracy.
        theta[["Z.var"]]<-ifelse(theta[["Z.var"]]<=.Machine$double.eps^0.5,
                                 sum(sweep(Z,2,apply(Z,2,mean))^2)/(n-d),
                                 theta[["Z.var"]])
        Z.pZK<-matrix(0,n,G)
        for(i in seq_len(n)) Z.pZK[i,cl[["cluster"]][i]]<-1
        
      }else{
        theta<-list(Z.mean = rbind(apply(Z,2,mean)),
                    Z.var = sum(sweep(Z,2,apply(Z,2,mean))^2)/(n-d),
                    Z.pK = 1)
      }
    }else theta<-resume
    
    E.step<-function(theta){
      Z.pZK<-with(theta,cbind(sapply(seq_len(G),function(g) Z.pK[g]*dmvnorm(Z,Z.mean[g,],Z.var[g]*diag(1,nrow=d)))))
      sweep(Z.pZK,1,apply(Z.pZK,1,sum),"/")
    }
    
    M.step<-function(Z.pZK){
      Z.pK <- apply(Z.pZK,2,mean)
      Z.mean <- sweep(crossprod(Z.pZK,Z),1,Z.pK,"/")/n
      ## Alternative implementation of the above vectorized implementation: may or may not be slower:
      ##t(sapply(seq_len(G),function(g) sapply(seq_len(d),function(j) weighted.mean(Z[,j],Z.pZK[,g]))))
      Z.var <- sapply(seq_len(G),function(g)
                      sum(Z.pZK[,g]*apply(sweep(Z,2,Z.mean[g,])^2,1,mean))/sum(Z.pZK[,g])
                      )
      list(Z.pK=Z.pK,Z.mean=Z.mean,Z.var=Z.var)
    }
    
    llk<-function(theta,Z.pZK){
      with(theta,
           sum(apply(Z,1,
                     function(z) log(sum(sapply(seq_len(G),
                                                function(g) Z.pK[g]*dmvnorm(z,Z.mean[g,],Z.var[g]*diag(1,nrow=d)))))
                     )
               )
           )
    }
    
    theta.old<-theta
    
    EMloop<-try(
                {
                  for(it in 1:EM.maxit){
                    converged<-FALSE
                    Z.pZK<-E.step(theta)
                    theta<-M.step(Z.pZK)
                    if(isTRUE(all.equal(theta.old,theta,tolerance=EM.tol))) {
                      converged=TRUE
                      break
                    }
                    theta.old<-theta
                  }   
                },
                silent=TRUE)

    Z.K<-apply(Z.pZK,1,which.max)
    
    if(inherits(EMloop,"try-error") || with(theta,max(Z.var)/min(Z.var))>.Machine$double.eps^-0.5){
      llk<-Inf
      bic<--Inf
    }else{
      llk<-llk(theta,Z.pZK)
      bic<--2*llk+(G-1 + d*G + G)*log(n)
      break
    }
  }
  return(with(theta,list(Z.mean=Z.mean,Z.var=Z.var,Z.K=Z.K,Z.pK=Z.pK,Z.pZK=Z.pZK,llk=llk,bic=bic,converged=converged,iterations=it)))
}
