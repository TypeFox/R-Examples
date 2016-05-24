interventionalInferenceAdvanced <-
function(y,X0,X1,cond=NULL,inhibition=NULL,inhibitors=NULL,max.indeg,g=NULL,Sigma=NULL,inferParents=NULL,allowSelfEdges=TRUE,perfect=FALSE,fixedEffect=FALSE,mechanismChange=FALSE,priorType="uninformed",priorGraph=NULL,priorStrength=3,fittedValues=FALSE) {
  # n is the number of samples.
  # P is the number of nodes.
  # y is an (n x P) matrix filled with response values.
  # X0 is an (n x a) matrix giving the part of the design matrix that is the same for all models. For no intercept, put X0=NULL.
  # X1 is an (n x P) matrix filled with the predictors.
  # Sigma is an (n x n) symmetric positive definite matrix giving the covariance of the responses (proportional to sigma^2).
  # cond is an (n x 1) matrix giving the condition of each sample, filled with the numbers 1,...,#conditions.
  # inhibition is a (#conditions x #inhibitors) binary matrix, where entry (c,i)=1 iff inhibitor i is active in condition c.
  # inhibitors is a (#inhibitors x P) matrix, where (i,p)=1 iff inhibitor i acts on node p (eg it removes edges OUT from node p).
  # inferParents is a vector of nodes for which to infer parents. If left blank, parents are inferred for all nodes. Any values in y for nodes that are not inferred are ignored.
  # allowSelfEdges if FALSE then self edges are forbidden, causing the diagonal of the pep matrix to be 0.
  # perfect: apply perfect interventions?
  # fixedEffect: apply fixed effect interventions?
  # mechanismChange: apply mechanism change interventions?
  # max.indeg: maximum indegree for each node.
  # g: the 'g' in Zellner's g-prior, by default set to be n.
  # priorType: type of prior to use. "uninformed", "Hamming" or "Mukherjee" ('only penalise unexpected' prior)
  # priorGraph: a (P x P) matrix specifying the prior network.
  # priorStrength: the prior strength parameter, ignored (but don't set it to NA) if priorGraph=NULL.
  # fittedValues: calculate fitted values?
  n<-dim(y)[1]
  P<-dim(y)[2]
  cat("n =",n,", nodes =",P,"\n")
  if (!is.null(X0) && dim(X0)[1]!=n) {stop("X0 must have dimension (n x a).\n")}
  if (dim(X1)[1]!=n | dim(X1)[2]!=P) {stop("X1 must have dimension (n x P).\n")}
  if (!is.null(Sigma) && dim(Sigma)[1]!=n && dim(Sigma)[2]!=n) {stop("Sigma must have dimension (n x n).\n")}
  if (perfect | fixedEffect | mechanismChange) {
    if (is.null(cond)) {stop("cond must be specified if an intervention model is used.\n")}
    if (is.null(inhibition)) {stop("inhibition must be specified if an intervention model is used.\n")}
    if (is.null(inhibitors)) {stop("inhibitors must be specified if an intervention model is used.\n")}
    cds<-dim(inhibition)[1] # number of conditions
    ins<-dim(inhibition)[2] # number of inhibitors
    if (length(cond)!=n) {stop("cond must have length n.\n")}
    if (max(cond)>cds) {stop("inhibition must have dimension (conditions x inhibitors).\n")}
    if (dim(inhibitors)[1]!=ins | dim(inhibitors)[2]!=P) {stop("inhibitors must have dimension (inhibitors x P).\n")}
  }
  if (is.null(inferParents)) {inferParents<-1:P}
  if (is.null(max.indeg) || max.indeg>P) {stop("max.indeg must be less than P.\n")}
  if (is.null(g)) {g<-n}
  if (priorType!="uninformed" & priorType!="Hamming" & priorType!="Mukherjee") {stop("priorType must be 'uninformed', 'Hamming' or 'Mukherjee'.\n")}
  if (priorType!="uninformed" & is.null(priorGraph)) {stop("priorGraph must be specified with an informed prior.\n")}
  if (!is.null(priorGraph) && (dim(priorGraph)[1]!=P | dim(priorGraph)[2]!=P)) {stop("priorGraph must have dimension (P x P).\n")}
  if (length(priorStrength)>1) {sort(priorStrength)}
  inputs<-list(y=y,X0=X0,X1=X1,cond=cond,inhibition=inhibition,inhibitors=inhibitors,max.indeg=max.indeg,g=g,Sigma=Sigma,inferParents=inferParents,
    allowSelfEdges=allowSelfEdges,perfect=perfect,fixedEffect=fixedEffect,mechanismChange=mechanismChange,
    priorType=priorType,priorGraph=priorGraph,priorStrength=priorStrength,n=n,P=P)
  # Part zero: prepare to remove covariance structure from everything.
  if (!is.null(Sigma)) {R<-t(chol(Sigma))}
  # Part one: remove component of X0 from y.
  if (!is.null(X0)) {
    a<-dim(X0)[2]
    cat("a =",a,"\n")
    #intercepts<-matrix(NA,a,P)
    if (is.null(Sigma)) {
      IP0<-diag(rep(1,n))-X0%*%solve(crossprod(X0),t(X0))
    } else {
      IP0<-solve(R)-solve(R,X0)%*%solve(crossprod(X0,solve(Sigma,X0)),t(solve(Sigma,X0)))
    }
    for (p in inferParents) {
      #intercepts[,p]<-solve(crossprod(X0),t(X0))%*%y[,p]
      y[,p]<-IP0%*%y[,p]
    }
  }
  # Part two: implement perfect interventions
  if (perfect) {
    for (j in 1:cds) {
      for (i in 1:ins) {
        if (inhibition[j,i]==1) {
          for (p in which(inhibitors[i,]==1)) {X1[which(cond==j),p]<-NA}
        }
      }
    }
  }
  # Part three: orthogonalise the predictors
  for (p in 1:P) {
    wh<-which(!is.na(X1[,p]))
    if (length(wh)==n) {
      X1[,p]<-IP0%*%X1[,p]
    } else {
      if (is.null(Sigma)) {
        X1[wh,p]<-(diag(rep(1,length(wh)))-X0[wh,]%*%solve(crossprod(X0[wh,]),t(X0[wh,])))%*%X1[wh,p]
      } else {
        X1[wh,p]<-(solve(R[wh,wh])-solve(R[wh,wh],X0[wh,])%*%solve(crossprod(X0[wh,],solve(Sigma[wh,wh],X0[wh,])),t(solve(Sigma[wh,wh],X0[wh,]))))%*%X1[wh,p]
      }
      X1[which(is.na(X1[,p])),p]<-0
    }
    if (sd(X1[,p])==0) {cat("Predictor",p,"is either constant or always inhibited.\n")}
  }
  # Part four: implement fixed effect interventions
  if (fixedEffect) {
    fe<-matrix(0,n,ins)
    for (i in 1:ins) {
      for (j in 1:cds) {
        if (inhibition[j,i]==1) {fe[which(cond==j),i]<-1}
      }

    }
    # Part four (a): orthogonalise fixed effects
    fe<-IP0%*%fe
  }
  # Part five: the prior
	grphs<-countGraphs(P,max.indeg)
	prior<-matrix(0,P,grphs)
  if (!is.null(priorGraph) & priorType!="uninformed") {
    if (priorType=="Hamming") {
      cat("Calculating Hamming prior distances...\n")
    } else if (priorType=="Mukherjee") {
      cat("Calculating Mukherjee prior distances...\n")
    } else {
      stop("Prior type not supported.\n")
    }
    parents<-rep(0,P)
    for (i in 1:countGraphs(P,max.indeg)) {
      wh<-which(parents==1)
      for (p in 1:P) {
        if (priorType=="Hamming") {# SHD prior
          prior[p,i]<-sum(abs(parents-priorGraph[,p]))
        } else if (priorType=="Mukherjee") {# OPU prior (Sach prior)
          prior[p,i]<-length(which(priorGraph[wh,p]==0))
        }
      }
      parents<-nxt(parents,max.indeg)
    }
  }
	# Part six: initialise
  ll<-matrix(NA,P,grphs)
  rownames(ll)<-colnames(X1)
  parentSets<-ll
	parents<-rep(0,P)
	st<-Sys.time()
	cat("Processing",grphs,"models",date(),"\n")
	# Part seven: The Main Loop
	for (m in 1:grphs) {
    parentSets[,m]<-parents # record parent set
    if (mechanismChange) {
      X<-matrix(NA,n,0)
      for (p in which(parents==1)) {
        unused.obs<-1:n
        for (i in c(which(inhibitors[,p]==1),0)) {
          if (i>0) {
            used.obs<-which(cond %in% which(inhibition[,i]==1))
          } else {
            used.obs<-unused.obs
          }
          if (length(used.obs)>length(intersect(used.obs,unused.obs))) {cat("Error: Mechanism change interventions do not support combinations of inhibitors that target the same protein.\n")}
          if (length(used.obs)>0) {
            newcol<-rep(0,n)
            newcol[used.obs]<-X1[used.obs,p]
            if (is.null(Sigma)) {
              newcol[used.obs]<-(diag(rep(1,length(used.obs)))-X0[used.obs,]%*%solve(crossprod(X0[used.obs,]),t(X0[used.obs,])))%*%newcol[used.obs]
            } else {
              newcol[used.obs]<-(solve(R[used.obs,used.obs])-solve(R[used.obs,used.obs],X0[used.obs,])%*%solve(crossprod(X0[used.obs,],solve(Sigma[used.obs,used.obs],X0[used.obs,])),t(solve(Sigma[used.obs,used.obs],X0[used.obs,]))))%*%newcol[used.obs]
            }
            if (sd(newcol)>0) {X<-cbind(X,newcol)} else {cat("sd zero\n")}
            unused.obs<-setdiff(unused.obs,used.obs)
          }
        }
      }
    } else {
      X<-matrix(X1[,which(parents==1)],n,sum(parents)) # put together design matrix for this model
    }
		if (fixedEffect) {
      for (i in 1:ins) {
        if (max(parents*inhibitors[i,])==1) {X<-cbind(X,fe[,i])}
      }
    }
    b<-dim(X)[2] # number of betas
    if (b==0) { # null model
      H<-diag(rep(1,n))
    } else {
      H<-diag(rep(1,n))-(g/(g+1))*X%*%solve(crossprod(X),t(X))
    }
    if (allowSelfEdges) {
  		for (p in inferParents) {
        ll[p,m]<--b/2*log(1+g)-(n-a)/2*log(crossprod(y[,p],H%*%y[,p]))
      }
    } else {
  		for (p in setdiff(inferParents,which(parents==1))) {
        ll[p,m]<--b/2*log(1+g)-(n-a)/2*log(crossprod(y[,p],H%*%y[,p]))
      }
    }    
		parents <- nxt(parents,max.indeg)
		if (m==100) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/100,"minutes.\n")}
		if (m==1000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/1000,"minutes.\n")}
		if (m==10000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/10000,"minutes.\n")}
		if (m==100000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/100000,"minutes.\n")}
	}
	cat("Actual duration",difftime(Sys.time(),st,units="mins"),"minutes.\n")
  # Part eight: renormalisation & MAP model
  cat("Renormalising...\n")
  parentCount<-apply(parentSets,2,sum)
  if (length(priorStrength)>1 & max(prior)>0) { # Perform Empirical Bayes to estimate the best one.
    marginal.likelihood<-matrix(0,P,length(priorStrength))
    ll.Inf<-ll
    ll.Inf[which(is.na(ll))]<--Inf
    st<-Sys.time()
 	  for (i in 1:length(priorStrength)) {
      # 'soft' multiplicity correction - prior is penalised and then renormalised
      normalisedPrior<-matrix(-log(1+max.indeg),P,grphs)# deal with the null model here.
      for (j in 1:max.indeg) {
        wh<-which(parentCount==j)
        normalisedPrior[,wh]<--log(length(wh))-log(1+max.indeg)
      }
      normalisedPrior <- normalisedPrior - priorStrength[i]*prior
      normalisedPrior <- normalisedPrior - log(apply(exp(normalisedPrior),1,sum))
      marginal.likelihood[,i] <- apply(ll+normalisedPrior,1,min,na.rm=TRUE)      
      marginal.likelihood[,i] <- marginal.likelihood[,i]+log(apply(exp(ll.Inf+normalisedPrior-marginal.likelihood[,i]),1,sum))
      if (i==20) {cat("Estimated duration of empirical Bayes calculation",difftime(Sys.time(),st,units="mins")*length(priorStrength)/20,"minutes.\n")}
    }
    log.ml.prod <- apply(marginal.likelihood,2,sum)
    if (length(which(log.ml.prod==max(log.ml.prod)))>1) {cat("Multiple solutions for Emperical Bayes - taking weakest prior strength.\n")}
    ebPriorStrength <- priorStrength[which.max(log.ml.prod)]
    const1 <- marginal.likelihood[,which.max(log.ml.prod)]
    # Now that we have chosen the prior strength, repeat the calculations.
    normalisedPrior<-matrix(-log(1+max.indeg),P,grphs)# deal with the null model here.
    for (j in 1:max.indeg) {
      wh<-which(parentCount==j)
      normalisedPrior[,wh]<--log(length(wh))-log(1+max.indeg)
    }
    normalisedPrior <- normalisedPrior-ebPriorStrength*prior
    normalisedPrior <- normalisedPrior-log(apply(exp(normalisedPrior),1,sum))
    lpost <- ll.Inf + normalisedPrior - const1
  } else { # No empirical Bayes
    # 'soft' multiplicity correction - prior is penalised and then renormalised
    normalisedPrior<-matrix(-log(1+max.indeg),P,grphs)# deal with the null model here.
    for (j in 1:max.indeg) {
      wh<-which(parentCount==j)
      normalisedPrior[,wh]<--log(length(wh))-log(1+max.indeg)
    }
    normalisedPrior<-normalisedPrior-priorStrength*prior
    normalisedPrior<-normalisedPrior-log(apply(exp(normalisedPrior),1,sum))
    const1<-apply(ll+normalisedPrior,1,min,na.rm=TRUE)
    ll[which(is.na(ll))]<--Inf
    const1<-const1+log(apply(exp(ll+normalisedPrior-const1),1,sum))
    lpost<-ll+normalisedPrior-const1
    marginal.likelihood <- const1
    ebPriorStrength <- NULL 
  }
  MAP<-matrix(NA,P,P)
  colnames(MAP)<-colnames(X1)
  rownames(MAP)<-colnames(X1)
  pep<-MAP
  MAPprob<-rep(NA,P)
  names(MAPprob)<-colnames(X1)
  MAPmodel<-MAPprob
  for (p in 1:P) {
    if (length(which.max(lpost[p,]))>0) {
      MAPmodel[p]<-which.max(lpost[p,])
      MAP[,p]<-parentSets[,MAPmodel[p]]
      MAPprob[p]<-exp(lpost[p,MAPmodel[p]])
    }
  }
  # Part nine: model averaging
  cat("Calculating posterior edge probabilities...\n")
	for (i in 1:P) {
    for (j in 1:P) {
      pep[i,j]<-sum(exp(lpost[j,which(parentSets[i,]==1)]))
		}
	}
  # Part ten: fitted values
  if (fittedValues) {
  	parents<-rep(0,P)
  	st<-Sys.time()
  	cat("Second pass to calculate fitted values.\n")
  	yhat<-matrix(0,n,P)
  	cat("Processing",grphs,"models",date(),"\n")
  	for (m in 1:grphs) {
      if (mechanismChange) {
        X<-matrix(NA,n,0)
        for (p in which(parents==1)) {
          unused.obs<-1:n
          for (i in c(which(inhibitors[,p]==1),0)) {
            if (i>0) {
              used.obs<-which(cond %in% which(inhibition[,i]==1))
            } else {
              used.obs<-unused.obs
            }
            if (length(used.obs)>length(intersect(used.obs,unused.obs))) {cat("Error: Mechanism change interventions do not support combinations of inhibitors that target the same protein.\n")}
            if (length(used.obs)>0) {
              newcol<-rep(0,n)
              newcol[used.obs]<-X1[used.obs,p]
              if (is.null(Sigma)) {
                newcol[used.obs]<-(diag(rep(1,length(used.obs)))-X0[used.obs,]%*%solve(crossprod(X0[used.obs,]),t(X0[used.obs,])))%*%newcol[used.obs]
              } else {
                newcol[used.obs]<-(solve(R[wh,wh])-solve(R[used.obs,used.obs],X0[used.obs,])%*%solve(crossprod(X0[used.obs,],solve(Sigma[used.obs,used.obs],X0[used.obs,])),t(solve(Sigma[used.obs,used.obs],X0[used.obs,]))))%*%newcol[used.obs]
              }
              if (sd(newcol)>0) {X<-cbind(X,newcol)} else {cat("sd zero\n")}
              unused.obs<-setdiff(unused.obs,used.obs)
            }
          }
        }
      } else {
        X<-matrix(X1[,which(parents==1)],n,sum(parents)) # put together design matrix for this model
      }
  		if (fixedEffect) {
        for (i in 1:ins) {
          if (max(parents*inhibitors[i,])==1) {X<-cbind(X,fe[,i])}
        }
      }
      b<-dim(X)[2] # number of betas
      if (b==0) { # null model
        H<-matrix(0,n,n)
      } else {
        H<-(g/(g+1))*X%*%solve(crossprod(X),t(X))
      }
  		for (p in inferParents) {
        yhat[,p]<-yhat[,p]+exp(lpost[p,m])*H%*%y[,p]
      }   
  		parents <- nxt(parents,max.indeg)
  		if (m==100) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/100,"minutes.\n")}
  		if (m==1000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/1000,"minutes.\n")}
  		if (m==10000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/10000,"minutes.\n")}
  		if (m==100000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/100000,"minutes.\n")}
  	}
  	cat("Actual duration",difftime(Sys.time(),st,units="mins"),"minutes.\n")
    # Map fitted values back onto the original scale
    if (is.null(Sigma)) {
      P0<-X0%*%solve(crossprod(X0),t(X0))
      for (p in inferParents) {
        yhat[,p]<-yhat[,p]+P0%*%inputs$y[,p]
      }
    } else {
      P0<-X0%*%solve(crossprod(X0,solve(Sigma,X0)),t(solve(Sigma,X0)))
      for (p in inferParents) {
        yhat[,p]<-R%*%yhat[,p]+P0%*%inputs$y[,p]
      }    
    }
  } else {
 	  yhat<-matrix(NA,n,P)
  }
  return(list(pep=pep,MAP=MAP,parentSets=parentSets,ll=ll,lpost=lpost,MAPprob=MAPprob,MAPmodel=MAPmodel,marginal.likelihood=marginal.likelihood,ebPriorStrength=ebPriorStrength,yhat=yhat,inputs=inputs))
}
