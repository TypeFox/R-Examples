interventionalInference <-
function(y,X0,X1,Z=NULL,max.indeg,g=NULL,Sigma=NULL,inferParents=NULL,allowSelfEdges=TRUE,perfectOut=FALSE,fixedEffectOut=FALSE,mechanismChangeOut=FALSE,perfectIn=FALSE,fixedEffectIn=FALSE,mechanismChangeIn=FALSE,priorType="uninformed",priorGraph=NULL,priorStrength=3,fittedValues=FALSE) {
  # n is the number of samples.
  # P is the number of nodes.
  # y is an (n x P) matrix filled with response values.
  # X0 is an (n x a) matrix giving the part of the design matrix that is the same for all models. For no intercept, put X0=NULL.
  # X1 is an (n x P) matrix filled with the predictors.
  # Z is an (n x P) binary matrix which is 1 iff node j is inhibited in sample i.
  # Sigma is an (n x n) symmetric positive definite matrix giving the covariance of the responses (proportional to sigma^2).
  # inferParents is a vector of nodes for which to infer parents. If left blank, parents are inferred for all nodes. Any values in y for nodes that are not inferred are ignored.
  # allowSelfEdges if FALSE then self edges are forbidden, causing the diagonal of the pep matrix to be 0.
  # perfectOut: apply perfect-out interventions?
  # fixedEffectOut: apply fixed-effect-out interventions?
  # mechanismChangeOut: apply mechanism-change-out interventions?
  # perfecIn: apply perfect-in interventions?
  # fixedEffectIn: apply fixed-effect-in interventions?
  # mechanismChangeIn: apply mechanism-change-in interventions?
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
  if (perfectOut | fixedEffectOut | mechanismChangeOut | perfectIn | fixedEffectIn | mechanismChangeIn) {
    if (is.null(Z)) {stop("Z must be specified if an intervention model is used.\n")}
  }
  if (is.null(inferParents)) {inferParents<-1:P}
  if (is.null(max.indeg) || max.indeg>P) {stop("max.indeg must be less than P.\n")}
  if (is.null(g)) {g<-n}
  if (priorType!="uninformed" & priorType!="Hamming" & priorType!="Mukherjee") {stop("priorType must be 'uninformed', 'Hamming' or 'Mukherjee'.\n")}
  if (priorType!="uninformed" & is.null(priorGraph)) {stop("priorGraph must be specified with an informed prior.\n")}
  if (!is.null(priorGraph) && (dim(priorGraph)[1]!=P | dim(priorGraph)[2]!=P)) {stop("priorGraph must have dimension (P x P).\n")}
  if (length(priorStrength)>1) {sort(priorStrength)}
  if ((perfectIn | perfectOut) & (mechanismChangeIn | mechanismChangeOut)) {stop("Perfect and mechanism change interventions cannot be used together.\n")}
  if (mechanismChangeIn & mechanismChangeOut) {stop("MC-in and MC-out not currently implemented.\n")}
  inputs<-list(y=y,X0=X0,X1=X1,Z=Z,max.indeg=max.indeg,g=g,Sigma=Sigma,inferParents=inferParents,
    allowSelfEdges=allowSelfEdges,perfectOut=perfectOut,fixedEffectOut=fixedEffectOut,mechanismChangeOut=mechanismChangeOut,
    perfectIn=perfectIn,fixedEffectIn=fixedEffectIn,mechanismChangeIn=mechanismChangeIn,
    priorType=priorType,priorGraph=priorGraph,priorStrength=priorStrength,n=n,P=P)
  # Part zero: prepare to remove covariance structure from everything.
  if (!is.null(Sigma)) {R<-t(chol(Sigma))}
  # Part one: remove component of X0 from y.
  if (!is.null(X0)) {
    a<-dim(X0)[2]
    cat("a =",a,"\n")
    for (p in inferParents) {
      if (perfectIn & max(Z[,p])==1) {
        obs<-which(Z[,p]==0)
      } else {
        obs<-1:n
      }
      X0p<-matrix(X0[obs,],length(obs),a)        
      if (is.null(Sigma)) {
        IP0<-diag(rep(1,length(obs)))-X0p%*%solve(crossprod(X0p),t(X0p))
      } else {
        IP0<-solve(R[obs,obs])-solve(R[obs,obs],X0p)%*%solve(crossprod(X0p,solve(Sigma[obs,obs],X0p)),t(solve(Sigma[obs,obs],X0p)))
      }
      y[obs,p]<-IP0%*%y[obs,p]
    }
  }
  # Part two: implement perfect-out interventions
  if (perfectOut) {X1[which(Z==1)]<-NA}
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
  # Part four: orthogonalise fixed-effect-out interventions
  if (fixedEffectOut) {fe<-IP0%*%Z}
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
    if (mechanismChangeOut) { # Mechanism-change-out interventions
      X<-matrix(NA,n,0)
      for (p in which(parents==1)) {
        if (max(Z[,p])==1) {
          newcols<-matrix(0,n,2)
          wh1<-which(Z[,p]==1)
          wh0<-which(Z[,p]==0)
          if (is.null(Sigma)) {
            newcols[wh1,1]<-(diag(rep(1,length(wh1)))-X0[wh1,]%*%solve(crossprod(X0[wh1,]),t(X0[wh1,])))%*%X1[wh1,p]
            newcols[wh0,2]<-(diag(rep(1,length(wh0)))-X0[wh0,]%*%solve(crossprod(X0[wh0,]),t(X0[wh0,])))%*%X1[wh0,p]
          } else {
            newcol[wh1,1]<-(solve(R[wh1,wh1])-solve(R[wh1,wh1],X0[wh1,])%*%solve(crossprod(X0[wh1,],solve(Sigma[wh1,wh1],X0[wh1,])),t(solve(Sigma[wh1,wh1],X0[wh1,]))))%*%X1[wh1,p] 
            newcol[wh0,2]<-(solve(R[wh0,wh0])-solve(R[wh0,wh0],X0[wh0,])%*%solve(crossprod(X0[wh0,],solve(Sigma[wh0,wh0],X0[wh0,])),t(solve(Sigma[wh0,wh0],X0[wh0,]))))%*%X1[wh0,p]
          }
          X<-cbind(X,newcols)
        } else {
          X<-cbind(X,X1[,p])
        }
      }
    } else {
      X<-matrix(X1[,which(parents==1)],n,sum(parents)) # put together design matrix for this model
    }
		if (fixedEffectOut) {X<-cbind(X,fe[,which(parents==1 & apply(Z,2,max)==1)])}
    b<-dim(X)[2] # number of betas
    if (b==0) { # null model
      H<-diag(rep(1,n))
    } else {
      H<-diag(rep(1,n))-(g/(g+1))*X%*%solve(crossprod(X),t(X))
    }
    inhibitedResponses<-NULL
    uninhibitedResponses<-inferParents
    if (!allowSelfEdges) {uninhibitedResponses<-intersect(uninhibitedResponses,which(parents==0))}
    if (perfectIn | fixedEffectIn | mechanismChangeIn) {
      inhibitedResponses<-intersect(uninhibitedResponses,which(apply(Z,2,max)==1))
      uninhibitedResponses<-setdiff(uninhibitedResponses,inhibitedResponses)
    }
    for (p in uninhibitedResponses) {
      ll[p,m]<--b/2*log(1+g)-(n-a)/2*log(crossprod(y[,p],H%*%y[,p]))
    }
    # Part 7a deal with -in interventions (and -out ones as well!)
    for (p in inhibitedResponses) {
      if (perfectIn) {
        obs<-which(Z[,p]==0)
      } else {
        obs<-1:n
      }
      if (mechanismChangeIn) {
        b<-length(which(parents==1))
        X<-matrix(0,n,2*b)
        wh1<-which(Z[,p]==1)
        wh0<-which(Z[,p]==0)
        if (is.null(Sigma)) {
          X0p<-matrix(X0[wh1,],length(wh1),a)
          X[wh1,1:b]<-(diag(rep(1,length(wh1)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X1[wh1,which(parents==1)]
          X0p<-matrix(X0[wh0,],length(wh0),a)
          X[wh0,(b+1):(2*b)]<-(diag(rep(1,length(wh0)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X1[wh0,which(parents==1)]
        } else {
          X0p<-matrix(X0[wh1,],length(wh1),a)
          X[wh1,1:b]<-(solve(R[wh1,wh1])-solve(R[wh1,wh1],X0p)%*%solve(crossprod(X0p,solve(Sigma[wh1,wh1],X0p)),t(solve(Sigma[wh1,wh1],X0p))))%*%X1[wh1,which(parents==1)]
          X0p<-matrix(X0[wh0,],length(wh0),a)
          X[wh0,(b+1):(2*b)]<-(solve(R[wh0,wh0])-solve(R[wh0,wh0],X0p)%*%solve(crossprod(X0p,solve(Sigma[wh0,wh0],X0p)),t(solve(Sigma[wh0,wh0],X0p))))%*%X1[wh0,which(parents==1)]
        }
      } else if (perfectOut) {
        X<-matrix(X1[obs,which(parents==1)],length(obs),sum(parents))
        for (pa in which(parents==1)) {
          if (max(Z[obs,pa])==1) {
            wh0<-obs[which(Z[obs,pa]==0)]
            X0p<-matrix(X0[wh0,],length(wh0),a)
            if (is.null(Sigma)) {
              X[wh0,which(which(parents==1)==pa)]<-(diag(rep(1,length(wh0)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X1[wh0,pa]
            } else {
              X[wh0,which(which(parents==1)==pa)]<-(solve(R[wh0,wh0])-solve(R[wh0,wh0],X0p)%*%solve(crossprod(X0p,solve(Sigma[wh0,wh0],X0p)),t(solve(Sigma[wh0,wh0],X0p))))%*%X1[wh0,pa]
            }            
            X[which(Z[obs,pa]==1),which(which(parents==1)==pa)]<-0
          }
        }            
      } else {
        X<-matrix(X1[obs,which(parents==1)],length(obs),sum(parents))
      }
      if (fixedEffectIn & fixedEffectOut) {
        X<-cbind(X,Z[obs,union(p,which(parents==1 & apply(Z[obs,],2,max)==1))])
      } else if (fixedEffectIn) {
        X<-cbind(X,Z[obs,p])
      } else if (fixedEffectOut) {
        X<-cbind(X,Z[obs,which(parents==1 & apply(Z[obs,],2,max)==1)])
      }         
      X0p<-matrix(X0[obs,],length(obs),a)
      if (is.null(Sigma)) {
        X<-(diag(rep(1,length(obs)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X
      } else {
        X<-(solve(R[obs,obs])-solve(R[obs,obs],X0p)%*%solve(crossprod(X0p,solve(Sigma[obs,obs],X0p)),t(solve(Sigma[obs,obs],X0p))))%*%X
      }
      b<-dim(X)[2] # number of betas
      if (b==0) { # null model
        H<-diag(rep(1,length(obs)))
      } else {
        H<-diag(rep(1,length(obs)))-(g/(g+1))*X%*%solve(crossprod(X),t(X))
      }
      ll[p,m]<--b/2*log(1+g)-(length(obs)-a)/2*log(crossprod(y[obs,p],H%*%y[obs,p]))      
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
      parentSets[,m]<-parents # record parent set
      if (mechanismChangeOut) { # Mechanism-change-out interventions
        X<-matrix(NA,n,0)
        for (p in which(parents==1)) {
          if (max(Z[,p])==1) {
            newcols<-matrix(0,n,2)
            wh1<-which(Z[,p]==1)
            wh0<-which(Z[,p]==0)
            if (is.null(Sigma)) {
              newcols[wh1,1]<-(diag(rep(1,length(wh1)))-X0[wh1,]%*%solve(crossprod(X0[wh1,]),t(X0[wh1,])))%*%X1[wh1,p]
              newcols[wh0,2]<-(diag(rep(1,length(wh0)))-X0[wh0,]%*%solve(crossprod(X0[wh0,]),t(X0[wh0,])))%*%X1[wh0,p]
            } else {
              newcol[wh1,1]<-(solve(R[wh1,wh1])-solve(R[wh1,wh1],X0[wh1,])%*%solve(crossprod(X0[wh1,],solve(Sigma[wh1,wh1],X0[wh1,])),t(solve(Sigma[wh1,wh1],X0[wh1,]))))%*%X1[wh1,p] 
              newcol[wh0,2]<-(solve(R[wh0,wh0])-solve(R[wh0,wh0],X0[wh0,])%*%solve(crossprod(X0[wh0,],solve(Sigma[wh0,wh0],X0[wh0,])),t(solve(Sigma[wh0,wh0],X0[wh0,]))))%*%X1[wh0,p]
            }
            X<-cbind(X,newcols)
          } else {
            X<-cbind(X,X1[,p])
          }
        }
      } else {
        X<-matrix(X1[,which(parents==1)],n,sum(parents)) # put together design matrix for this model
      }
  		if (fixedEffectOut) {X<-cbind(X,fe[,which(parents==1 & apply(Z,2,max)==1)])}
      b<-dim(X)[2] # number of betas
      if (b==0) { # null model
        H<-matrix(0,n,n)
      } else {
        H<-(g/(g+1))*X%*%solve(crossprod(X),t(X))
      }
      inhibitedResponses<-NULL
      uninhibitedResponses<-inferParents
      if (!allowSelfEdges) {uninhibitedResponses<-intersect(uninhibitedResponses,which(parents==0))}
      if (perfectIn | fixedEffectIn | mechanismChangeIn) {
        inhibitedResponses<-intersect(uninhibitedResponses,which(apply(Z,2,max)==1))
        uninhibitedResponses<-setdiff(uninhibitedResponses,inhibitedResponses)
      }
      for (p in uninhibitedResponses) {
        yhat[,p]<-yhat[,p]+exp(lpost[p,m])*H%*%y[,p]
      }
      # Part 10a deal with -in interventions (and -out ones as well!)
      for (p in inhibitedResponses) {
        if (perfectIn) {
          obs<-which(Z[,p]==0)
        } else {
          obs<-1:n
        }
        if (mechanismChangeIn) {
          b<-length(which(parents==1))
          X<-matrix(0,n,2*b)
          wh1<-which(Z[,p]==1)
          wh0<-which(Z[,p]==0)
          if (is.null(Sigma)) {
            X0p<-matrix(X0[wh1,],length(wh1),a)
            X[wh1,1:b]<-(diag(rep(1,length(wh1)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X1[wh1,which(parents==1)]
            X0p<-matrix(X0[wh0,],length(wh0),a)
            X[wh0,(b+1):(2*b)]<-(diag(rep(1,length(wh0)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X1[wh0,which(parents==1)]
          } else {
            X0p<-matrix(X0[wh1,],length(wh1),a)
            X[wh1,1:b]<-(solve(R[wh1,wh1])-solve(R[wh1,wh1],X0p)%*%solve(crossprod(X0p,solve(Sigma[wh1,wh1],X0p)),t(solve(Sigma[wh1,wh1],X0p))))%*%X1[wh1,which(parents==1)]
            X0p<-matrix(X0[wh0,],length(wh0),a)
            X[wh0,(b+1):(2*b)]<-(solve(R[wh0,wh0])-solve(R[wh0,wh0],X0p)%*%solve(crossprod(X0p,solve(Sigma[wh0,wh0],X0p)),t(solve(Sigma[wh0,wh0],X0p))))%*%X1[wh0,which(parents==1)]
          }
        } else if (perfectOut) {
          X<-matrix(X1[obs,which(parents==1)],length(obs),sum(parents))
          for (pa in which(parents==1)) {
            if (max(Z[obs,pa])==1) {
              wh0<-obs[which(Z[obs,pa]==0)]
              X0p<-matrix(X0[wh0,],length(wh0),a)
              if (is.null(Sigma)) {
                X[wh0,which(which(parents==1)==pa)]<-(diag(rep(1,length(wh0)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X1[wh0,pa]
              } else {
                X[wh0,which(which(parents==1)==pa)]<-(solve(R[wh0,wh0])-solve(R[wh0,wh0],X0p)%*%solve(crossprod(X0p,solve(Sigma[wh0,wh0],X0p)),t(solve(Sigma[wh0,wh0],X0p))))%*%X1[wh0,pa]
              }            
              X[which(Z[obs,pa]==1),which(which(parents==1)==pa)]<-0
            }
          }            
        } else {
          X<-matrix(X1[obs,which(parents==1)],length(obs),sum(parents))
        }
        if (fixedEffectIn & fixedEffectOut) {
          X<-cbind(X,Z[obs,union(p,which(parents==1 & apply(Z[obs,],2,max)==1))])
        } else if (fixedEffectIn) {
          X<-cbind(X,Z[obs,p])
        } else if (fixedEffectOut) {
          X<-cbind(X,Z[obs,which(parents==1 & apply(Z[obs,],2,max)==1)])
        }         
        X0p<-matrix(X0[obs,],length(obs),a)
        if (is.null(Sigma)) {
          X<-(diag(rep(1,length(obs)))-X0p%*%solve(crossprod(X0p),t(X0p)))%*%X
        } else {
          X<-(solve(R[obs,obs])-solve(R[obs,obs],X0p)%*%solve(crossprod(X0p,solve(Sigma[obs,obs],X0p)),t(solve(Sigma[obs,obs],X0p))))%*%X
        }
        b<-dim(X)[2] # number of betas
        if (b==0) { # null model
          H<-matrix(0,length(obs),length(obs))
        } else {
          H<-(g/(g+1))*X%*%solve(crossprod(X),t(X))
        }
        yhat[obs,p]<-yhat[obs,p]+exp(lpost[p,m])*H%*%y[obs,p]      
      }       
  		parents <- nxt(parents,max.indeg)
  		if (m==100) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/100,"minutes.\n")}
  		if (m==1000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/1000,"minutes.\n")}
  		if (m==10000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/10000,"minutes.\n")}
  		if (m==100000) {cat("Estimated duration",difftime(Sys.time(),st,units="mins")*grphs/100000,"minutes.\n")}
  	}
  	cat("Actual duration",difftime(Sys.time(),st,units="mins"),"minutes.\n")
    # Map fitted values back onto the original scale
    for (p in inferParents) {
      if (perfectIn & max(Z[,p])==1) {
        obs<-which(Z[,p]==0)
        yhat[which(Z[,p]==1),p]<-NA
      } else {
        obs<-1:n
      }
      X0p<-matrix(X0[obs,],length(obs),a)        
      if (is.null(Sigma)) {
        P0<-diag(rep(1,length(obs)))-X0p%*%solve(crossprod(X0p),t(X0p))
        yhat[obs,p]<-yhat[obs,p]+P0%*%inputs$y[obs,p]
      } else {
        P0<-solve(R[obs,obs],X0p)%*%solve(crossprod(X0p,solve(Sigma[obs,obs],X0p)),t(solve(Sigma[obs,obs],X0p)))
        yhat[obs,p]<-R[obs,obs]%*%yhat[obs,p]+P0%*%inputs$y[obs,p]
      }
    }
  } else {
 	  yhat<-matrix(NA,n,P)
  }
  return(list(pep=pep,MAP=MAP,parentSets=parentSets,ll=ll,lpost=lpost,MAPprob=MAPprob,MAPmodel=MAPmodel,marginal.likelihood=marginal.likelihood,ebPriorStrength=ebPriorStrength,yhat=yhat,inputs=inputs))
}
