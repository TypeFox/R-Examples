oglmx<-function(formulaMEAN,formulaSD=NULL,data,start=NULL,weights=NULL,link="probit",constantMEAN=TRUE,constantSD=TRUE,beta=NULL,delta=NULL,threshparam=NULL,analhessian=TRUE,sdmodel=expression(exp(z)),SameModelMEANSD=FALSE,na.action=TRUE,savemodelframe=FALSE,Force=FALSE,robust=FALSE){  
  call<-match.call()
  # formulaMEAN specifies a formula for the conditional mean of the latent variable
  # formulaSD specifies a formula for the standard deviation of the error term of the latent variable. If NULL a homoskedastic
  # model is estimated. This may be a standard ordered probit or interval regression  since cutoffs can be specified and the constant in the equation for the mean included.
  # inputs of beta, delta and threshparam allow the user to prespecify the values of particular elements of the parameter. I.e. the standard ordered probit sets beta[1]==0 and delta[1]==0 if sdmodel=expression(exp(z))  
  # if specified beta should be of length equal to the number of parameters in beta, delta the number of parameters in delta, threshparam the number of outcomes - 1. Write NA for a parameter to be estimated.
  # constantMEAN and constantSD indicate that the equations for the mean of the latent variable and the standard deviation of the error include a constant.

  
# Extract necessary data
dataframeMEAN<-model.frame(formulaMEAN,data=data,na.action=na.pass)
termsMEAN<-terms(dataframeMEAN)
X<-model.matrix(formulaMEAN,dataframeMEAN)
Y<-model.response(dataframeMEAN)
KeepY<-!is.na(Y)
KeepX<-!apply(X,1,.IsNARow)
if (!is.null(formulaSD) & !SameModelMEANSD){
  dataframeSD<-model.frame(formulaSD,data=data,na.action=na.pass)
  Z<-model.matrix(formulaSD,dataframeSD)
  KeepZ<-!apply(Z,1,.IsNARow)
  X<-X[KeepX & KeepZ & KeepY, ,drop=FALSE]
  Z<-Z[KeepX & KeepZ & KeepY, ,drop=FALSE]
  Y<-Y[KeepX & KeepZ & KeepY]
  termsSD<-terms(dataframeSD)
} else {
  KeepZ<-!vector("logical",length(X))
  X<-X[KeepX & KeepY, ,drop=FALSE]
  Y<-Y[KeepX & KeepY]
  termsSD<-NULL
}
if (!is.null(weights)){
  weights<-weights[KeepX & KeepY & KeepZ]
  robust<-TRUE
  X<-X[!is.na(weights), ,drop=FALSE]
  Y<-Y[!is.na(weights)]
  if (!is.null(formulaSD) & !SameModelMEANSD){
    Z<-Z[!is.na(weights)]
  }
  weights<-weights[!is.na(weights)]
  weighted<-TRUE
  weights<-weights*length(weights)/sum(weights)
} else {
  weights<-rep_len(1,length(Y))
  weighted<-FALSE
}

NoVarModData<-data.frame(Y,weights)
  
No.Obs<-length(Y)
# calculate log-likelihood if probabilities are derived from proportions in sample
# used for calculation of McFadden's R^2
# BaselineLL<-.BaseLL(Y,weights)
termsMODEL<-list(termsMEAN,termsSD)
formulaMODEL<-list(formulaMEAN,formulaSD)
# If specified that there is no constant in the model remove it from the data frame.
if (!constantMEAN){
  savecolnames<-colnames(X)[colnames(X)!="(Intercept)"]
  X<-X[,colnames(X)!="(Intercept)",drop=FALSE]
  colnames(X)<-savecolnames
}

# collect variable means and check which variables are binary.
XVarMeans<-apply(X,2,mean)
XVarBinary<-apply(X,2,.checkbinary)  

  Heteroskedastic<-TRUE # set model to heteroskedastic, check call and then switch if not.
  if (!is.null(formulaSD)){
    # will check if the formula for the mean is identical to the 
    # the right hand side variables of a formula are accessible via terms in the attribute term.labels
    meaneqnames<-attr(terms(formulaMEAN),"term.labels")
    sdeqnames<-attr(terms(formulaSD),"term.labels")
    if (sum(is.na(match(meaneqnames,sdeqnames)))==sum(is.na(match(sdeqnames,meaneqnames))) & sum(is.na(match(sdeqnames,meaneqnames)))==0){
      if (constantSD==constantMEAN){
        SameModelMEANSD<-TRUE
        # collect the names and column numbers of variables that are in both the mean and variance equation
        meanandvarNAME<-colnames(X)
        meanandvarLOC<-c(1:ncol(X))
        meanandvarLOCZ<-meanandvarLOC<-meanandvarLOC[meanandvarNAME!="(Intercept)"]
        meanandvarNAME<-meanandvarNAME[meanandvarNAME!="(Intercept)"]
        BothMeanVar<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors=FALSE)
        ZVarMeans<-XVarMeans
        ZVarBinary<-XVarBinary
      } 
    }
  }

  if (!is.null(formulaSD) & !SameModelMEANSD){   
    if (!constantSD){Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]}
    ZVarMeans<-apply(Z,2,mean)
    ZVarBinary<-apply(Z,2,.checkbinary)
    # collect the names and column numbers of variables that are in both the mean and variance equation
    # find the colnames of Z that are the same as the colnames of X
    meanandvarLOC<-c(1:ncol(X))[!is.na(match(colnames(X),colnames(Z)))]
    # find the columns of Z that are in Z and X
    meanandvarLOCZ<-match(colnames(X),colnames(Z))[!is.na(match(colnames(X),colnames(Z)))]
    meanandvarNAME<-colnames(X)[meanandvarLOC]
    meanandvarLOC<-meanandvarLOC[meanandvarNAME!="(Intercept)"]
    meanandvarLOCZ<-meanandvarLOCZ[meanandvarNAME!="(Intercept)"]
    meanandvarNAME<-meanandvarNAME[meanandvarNAME!="(Intercept)"]
    BothMeanVar<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors=FALSE)
  } else if (is.null(formulaSD) & !SameModelMEANSD){
    Z<-as.matrix(rep(1,nrow(X)),ncol=1)
    ZVarMeans<-NULL
    ZVarBinary<-NULL
    Heteroskedastic<-FALSE
    if (is.null(delta) & is.null(threshparam)){
      # if no formula for the standard deviation is given, the threshold parameters are not specified or the standard deviation is not specified then use the unit variance assumption.
      calcdelta<-function(x){eval({z<-x;sdmodel})-1} 
      delta<-uniroot(calcdelta,c(-10,10),extendInt="yes",tol=.Machine$double.eps)$root # solve for delta to get the unit variance assumption
    }
    BothMeanVar=NULL
  }
  
  no.Xvar<-ncol(X)
  if (!SameModelMEANSD){no.Zvar<-ncol(Z)} else {no.Zvar<-no.Xvar}
  # count number of different outcomes and store.
  listoutcomes<-as.numeric(levels(as.factor(Y)))[order(as.numeric(levels(as.factor(Y))))]
  no.outcomes<-length(listoutcomes)
  if (no.outcomes>20 & Force==FALSE){
    stop("More than 20 different values for outcome variable.\n If you are sure you wish to estimate this model rerun command with Force option set to TRUE.")
  }
  # set the prespecified and not prespecified parts.
  # beta
  if (!is.null(beta) & length(beta)==1){
    storebeta<-beta
    beta<-rep(NA,no.Xvar)
    beta[1]<-storebeta
  } else if (!is.null(beta) & length(beta)>1){
    # check that the specified vector is of correct length
    if (length(beta)!=no.Xvar){stop("Specified beta vector of incorrect length.")}
  } else if (is.null(beta)){
    beta<-rep(NA,no.Xvar)
  }
  # delta
  if (!is.null(delta) & length(delta)==1){
    storedelta<-delta
    delta<-rep(NA,no.Zvar)
    delta[1]<-storedelta
  } else if (!is.null(delta) & length(delta)>1){
    # check that the specified vector is of correct length
    if (length(delta)!=no.Zvar){stop("Specified delta vector of incorrect length.")}
  } else if (is.null(delta)){
    delta<-rep(NA,no.Zvar)
  }
  # threshparam
  if (!is.null(threshparam) & length(threshparam)==1){
    storethreshparam<-threshparam
    threshparam<-rep(NA,no.outcomes-1)
    threshparam[1]<-storethreshparam
  } else if (!is.null(threshparam) & length(threshparam)>1){
    # check that the specified vector is of correct length
    if (length(threshparam)!=no.outcomes-1){stop("Specified vector of threshold parameters of incorrect length.")}
  } else if (is.null(threshparam)){
    threshparam<-rep(NA,no.outcomes-1)
  }

if (!is.null(beta)){
  collectbeta<-is.na(beta)
} else {
  collectbeta<-!vector("logical",no.Xvar)
}
if (!is.null(delta)){
  collectdelta<-is.na(delta)
} else {
    collectdelta<-!vector("logical",no.Zvar)
}
if (!is.null(threshparam)){
  collectthreshparam<-is.na(threshparam)
} else {
  collectthreshparam<-!vector("logical",no.outcomes-1)
}

no.betaparams<-sum(collectbeta)
no.deltaparams<-sum(collectdelta)
no.threshparams<-sum(collectthreshparam)
no.parameters<-no.betaparams+no.deltaparams+no.threshparams
Est.Parameters<-list(beta=collectbeta,delta=collectdelta,alpha=collectthreshparam)

  # specify the start vector in the case that it is given as null, give error message if not of correct length.
  if (is.null(start)){
    # start with a vector of zeros for the betas
    start<-vector("numeric",no.parameters)
    # if none of the delta parameters are set, set the first element so that the initial standard deviation is 0.5
    calcstartdelta<-function(x){eval({z<-x;sdmodel})-0.5}
    startdelta<-uniroot(calcstartdelta,c(-10,10),extendInt="yes")$root
    if (no.deltaparams==length(delta)){start[no.betaparams+1]<-startdelta}
    # more complicated for threshparam, should respect the order of the prespecified values
    if (no.threshparams>0){
      cutoff<-1
      for (i in 1:length(threshparam)){
        if (collectthreshparam[i]){
          start[no.betaparams+no.deltaparams+cutoff]<-(listoutcomes[i]+listoutcomes[i+1])/2
          cutoff<-cutoff+1
        } 
      }
    }
  } else if (length(start)!=no.parameters){
    stop("Specified vector of start values for parameters of incorrect length.")
  }
  # create the vectors of names matched to parameter estimates
  collectXnames<-colnames(X)[collectbeta]
  if (!SameModelMEANSD){
    collectZnames<-colnames(Z)[collectdelta]
  } else {
    collectZnames<-colnames(X)[collectdelta]
  }
  
  outputnames<-vector("character",0)
  if (no.betaparams>0){
    outputnames<-c(outputnames,collectXnames)
  }
  if (no.deltaparams>0){
    outputnames<-c(outputnames,collectZnames)
  }
  if (no.threshparams>0){
    collectnumbers<-c(1:length(threshparam))[collectthreshparam]
    outputnames<-c(outputnames,sapply(collectnumbers,function(x){paste("Threshold (", listoutcomes[x],"->",listoutcomes[x+1],")",sep="")}))
  }
  
  # save model frame if requested
  if (savemodelframe){
    modelframes<-list(X)
    if (!SameModelMEANSD){
      modelframes[[2]]<-Z
    } else {
      modelframes[[2]]<-X
    }
  } else {
    modelframes<-NULL
  }
  
  # break X (and Z) into chunks, relevant for each outcome. Delete X (and Z) and work with the already separated chunks
  # will reduce memory overhead in loop
  Xr<-list()
  if (!SameModelMEANSD){
    Zr<-list()
  } else {
    Zr<-NULL
  }

  Xr<-split.data.frame(X,Y,drop=FALSE)
  weightsr<-split(weights,Y,drop=FALSE)
  if (!SameModelMEANSD){
    Zr<-split.data.frame(Z,Y,drop=FALSE)
    rm(Z)
  }
  rm(X,Y)
  
# define function used to calculate likelihood, etc.
if (link=="logit"){
  ProbFunc<-function(p){plogis(p)}
  ProbFuncD<-function(p){dlogis(p)}
  ProbFuncDD<-function(p){dlogis(p)*(1-2*plogis(p))}
}
if (link=="probit"){
  ProbFunc<-function(p){pnorm(p)}
  ProbFuncD<-function(p){dnorm(p)}
  ProbFuncDD<-function(p){-p*dnorm(p)}
}
if (link=="cauchit"){
  ProbFunc<-function(p){pcauchy(p)}
  ProbFuncD<-function(p){dcauchy(p)}
  ProbFuncDD<-function(p){-2*p*dcauchy(p)/(1+p^2)}
}
if (link=="loglog"){
  ProbFunc<-function(p){exp(-exp(-p))}
  ProbFuncD<-function(p){exp(-(exp(-p)+p))}
  ProbFuncDD<-function(p){exp(-(exp(-p)+p))*(1+exp(-p))}
}
if (link=="cloglog"){
  ProbFunc<-function(p){1-exp(-exp(p))}
  ProbFuncD<-function(p){exp(p-exp(p))}
  ProbFuncDD<-function(p){exp(p-exp(p))*(1-exp(p))}
}
  
# function that calculates log likelihood, gradient and hessian given a set of parameter values
LLoglmx<-function(param,beta=NULL,delta=NULL,threshparam=NULL,analhessian=FALSE,robustmatrix=FALSE){
  if (is.null(beta)){ # if elements of the vector beta are not prespecified then all are included in param
    beta<-param[1:no.Xvar]
    param<-param[(no.Xvar+1):length(param)] # remove from param the elements allocated to the beta vector
  } else { # if not then fill NAs in beta with the first elements in param
    countNA<-sum(is.na(beta))
    if (countNA>0){beta[is.na(beta)]<-param[1:countNA]}
    param<-param[(countNA+1):length(param)]
  }
  # repeat to fill the delta vector
  if (is.null(delta)){ # if elements of the vector delta are not prespecified then all are included in param
    if (SameModelMEANSD){
      delta<-param[1:no.Xvar]
      param<-param[(no.Xvar+1):length(param)] # remove from param the elements allocated to the beta vector
    } else {
      delta<-param[1:no.Zvar]
      param<-param[(no.Zvar+1):length(param)] # remove from param the elements allocated to the beta vector
    }
  } else { # if not then fill NAs in delta with the first elements in param
    countNA<-sum(is.na(delta))
    if (countNA>0){delta[is.na(delta)]<-param[1:countNA]}
    param<-param[(countNA+1):length(param)]
  }
  # repeat to fill threshparam vector.
  if (is.null(threshparam)){
    threshparam<-param
  } else {
    if (length(param)>0){threshparam[is.na(threshparam)]<-param} else {stop("Insufficient number of parameters specified.")}
  }
  threshparam<-c(-Inf,threshparam,Inf)
  
  calcprobs<-function(outcome){
    # function that calculates relevant probabilities relevant to outcome
    # only to be called inside LLoglmx
    Xb<-Xr[[outcome]]%*%beta
    if (!SameModelMEANSD){
      Zdinv<-1/eval({z<-Zr[[outcome]]%*%delta;sdmodel})
    } else {
      Zdinv<-1/eval({z<-Xr[[outcome]]%*%delta;sdmodel})
    }
    Probs<-ProbFunc((threshparam[outcome+1]-Xb)*Zdinv)-ProbFunc((threshparam[outcome]-Xb)*Zdinv)
  }
  
  sdmodfirstderiv<-D(sdmodel,"z")
  sdmodsecondderiv<-D(sdmodfirstderiv,"z")
  
  delta<-as.matrix(delta)
  beta<-as.matrix(beta)
  
  vectorsprobs<-lapply(c(1:no.outcomes),calcprobs)
  # if weights are used then the standard log likelihood is no longer a relevant measure of model suitability.
  # need to calculate a pseudo-log likelihood, also for the baseline log-likelihood should take account of weights
  
  #loglikelihood<-sum(sapply(suppressWarnings(lapply(vectorsprobs,log)),sum))
  wloglikelihoodvecs<-list()
  for (i in 1:no.outcomes){
    wloglikelihoodvecs[[i]]<-suppressWarnings(log(vectorsprobs[[i]]))*weightsr[[i]]
  }
  loglikelihood<-sum(sapply(wloglikelihoodvecs,sum))
  
  # write function to work with the matrices for each outcome separately
  # produce the relevant gradient and hessian.
  # afterwards can sum. Allows the use of the Map function.
  
  getLLgradhess<-function(X,Z,w,index){
    j<-index
    Xb<-X%*%beta
    Zd<-Z%*%delta
    Zdinv<-1/eval({z<-Zd;sdmodel})
    if (j==1){
      frac0<-rep(-Inf,length(Zd))
      frac1<-(threshparam[j+1]-Xb)*Zdinv
    } else if (j==no.outcomes){
      frac1<-rep(Inf,length(Zd))
      frac0<-(threshparam[j]-Xb)*Zdinv
    } else {
      frac1<-(threshparam[j+1]-Xb)*Zdinv
      frac0<-(threshparam[j]-Xb)*Zdinv
    }
    
    # functions used to calculate score and hessian.
    calcscorebeta<-function(beta){
      # only to be called inside LLoglmx
      if (j==1){
        ProbderivBeta<-X[,beta]*Zdinv*(-ProbFuncD(frac1))/vectorsprobs[[j]]
      } else if (j==no.outcomes){
        ProbderivBeta<-X[,beta]*Zdinv*(ProbFuncD(frac0))/vectorsprobs[[j]]
      } else {
        ProbderivBeta<-X[,beta]*Zdinv*(ProbFuncD(frac0)-ProbFuncD(frac1))/vectorsprobs[[j]]
      }
      ProbderivBeta
    }
    
    calcscoredelta<-function(delta){
      # only to be called inside LLoglmx
      if (j==1){
        ProbderivDelta<- -Z[,delta]*eval({z<-Zd;sdmodfirstderiv})*Zdinv*frac1*ProbFuncD(frac1)/vectorsprobs[[j]]
      } else if (j==no.outcomes){
        ProbderivDelta<- Z[,delta]*eval({z<-Zd;sdmodfirstderiv})*Zdinv*frac0*ProbFuncD(frac0)/vectorsprobs[[j]]
      } else {
        ProbderivDelta<- -Z[,delta]*eval({z<-Zd;sdmodfirstderiv})*Zdinv*(frac1*ProbFuncD(frac1)-frac0*ProbFuncD(frac0))/vectorsprobs[[j]]
      }
      ProbderivDelta
    }
    
    calcscorethreshparam<-function(alpha){
      # only to be called inside LLoglmx
      if (alpha==j){
        Probderivalpha<- -Zdinv*ProbFuncD(frac0)/vectorsprobs[[j]] 
      } else if (alpha==j+1){
        Probderivalpha<- Zdinv*ProbFuncD(frac1)/vectorsprobs[[j]]
      } else {
        Probderivalpha<-vector("numeric",nrow(X))
      }
      Probderivalpha
    }
    
    scorevector<-vector("numeric",no.parameters)
    
    if (no.betaparams>0){
      vectorprobderivbeta<-sapply(c(1:length(beta))[collectbeta],calcscorebeta)
      scorevector[1:sum(collectbeta)]<-scorevector[1:sum(collectbeta)]+apply(vectorprobderivbeta*w,2,sum)
    }
    
    if (no.deltaparams>0){
      vectorprobderivdelta<-sapply(c(1:length(delta))[collectdelta],calcscoredelta)
      scorevector[(1+no.betaparams):(no.betaparams+no.deltaparams)]<-scorevector[(1+no.betaparams):(no.betaparams+no.deltaparams)]+apply(vectorprobderivdelta*w,2,sum)  
    }
    
    if (no.threshparams>0){
      vectorprobderivthreshparam<-sapply(c(2:(length(threshparam)-1))[collectthreshparam],calcscorethreshparam)
      scorevector[(1+no.betaparams+no.deltaparams):no.parameters]<-scorevector[(1+no.betaparams+no.deltaparams):no.parameters]+apply(vectorprobderivthreshparam*w,2,sum) 
    }
    
    if (robustmatrix){
      if (no.betaparams>0 & no.deltaparams>0){
        scorevecs<-cbind(vectorprobderivbeta,vectorprobderivdelta)
        if (no.threshparams>0){
          scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
        }
      } else if (no.betaparams>0){
        scorevecs<-vectorprobderivbeta
        if (no.threshparams>0){
          scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
        }
      } else if (no.deltaparams>0){
        scorevecs<-vectorprobderivdelta
        if (no.threshparams>0){
          scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
        }
      }
      BHHHmatrix<-matrix(vector("numeric",ncol(scorevecs)^2),nrow=ncol(scorevecs),ncol=ncol(scorevecs))
      for (v in 1:nrow(scorevecs)){
        obmat<-scorevecs[v,]%*%t(scorevecs[v,])
        BHHHmatrix<-BHHHmatrix+(w[v]^2)*obmat
      }
    } else {
      BHHHmatrix<-NULL
    }
    
    
    if (analhessian){
        calc2ndderivprobbetabeta<-function(x){ # x is a two element vector specifying location of each coefficient
          # only to be called inside LLoglmx
          if (j==1){
            probderiv2beta2<-sum((X[,x[1]]*X[,x[2]]*(Zdinv^2)*(ProbFuncDD(frac1))/vectorsprobs[[j]])*w)
          } else if (j==no.outcomes){
            probderiv2beta2<-sum((X[,x[1]]*X[,x[2]]*(Zdinv^2)*(-ProbFuncDD(frac0))/vectorsprobs[[j]])*w)
          } else {
            probderiv2beta2<-sum((X[,x[1]]*X[,x[2]]*(Zdinv^2)*(ProbFuncDD(frac1)-ProbFuncDD(frac0))/vectorsprobs[[j]])*w)
          }
          probderiv2beta2
        }
        
        calc2ndderivprobbetadelta<-function(x){ # x is a two element vector specifying location of each coefficient, beta 1st, delta 2nd
          # only to be called inside LLoglmx
          if (j==1){
            probderiv2betadelta<- X[,x[1]]*Z[,x[2]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac1)+frac1*ProbFuncDD(frac1))/vectorsprobs[[j]]
          } else if (j==no.outcomes){
            probderiv2betadelta<- X[,x[1]]*Z[,x[2]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(-ProbFuncD(frac0)-frac0*ProbFuncDD(frac0))/vectorsprobs[[j]]
          } else {
            probderiv2betadelta<- X[,x[1]]*Z[,x[2]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac1)-ProbFuncD(frac0)+frac1*ProbFuncDD(frac1)-frac0*ProbFuncDD(frac0))/vectorsprobs[[j]]
          }
          sum(probderiv2betadelta*w)
        }
        
        calc2ndderivprobbetaalpha<-function(x){ # x is a two element vector specifying location of each coefficient, beta 1st, alpha 2nd
          # only to be called inside LLoglmx
          if ((x[2]==j) & (j>1)){
            probderiv2betaalpha<- X[,x[1]]*(Zdinv^2)*ProbFuncDD(frac0)/vectorsprobs[[j]]
          } else if ((x[2]==j+1) & (j<no.outcomes)){
            probderiv2betaalpha<- -X[,x[1]]*(Zdinv^2)*ProbFuncDD(frac1)/vectorsprobs[[j]]
          } else {
            probderiv2betaalpha<-vector("numeric",nrow(Xr[[j]]))
          }
          sum(probderiv2betaalpha*w)
        }
        
        calc2ndderivprobdeltadelta<-function(x){ # x is a two element vector specifying location of each coefficient
          # only to be called inside LLoglmx
          if (j==1){
            probderiv2deltadelta<- Z[,x[1]]*Z[,x[2]]*Zdinv*((2*(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv-eval({z<-Zd;sdmodsecondderiv}))*(frac1*ProbFuncD(frac1))+(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv*((frac1^2)*ProbFuncDD(frac1)))/vectorsprobs[[j]]
          } else if (j==no.outcomes){
            probderiv2deltadelta<- Z[,x[1]]*Z[,x[2]]*Zdinv*((2*(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv-eval({z<-Zd;sdmodsecondderiv}))*(-frac0*ProbFuncD(frac0))+(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv*(-(frac0^2)*ProbFuncDD(frac0)))/vectorsprobs[[j]]
          } else {
            probderiv2deltadelta<- Z[,x[1]]*Z[,x[2]]*Zdinv*((2*(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv-eval({z<-Zd;sdmodsecondderiv}))*(frac1*ProbFuncD(frac1)-frac0*ProbFuncD(frac0))+(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv*((frac1^2)*ProbFuncDD(frac1)-(frac0^2)*ProbFuncDD(frac0)))/vectorsprobs[[j]]
          }
          sum(probderiv2deltadelta*w)
        }
        
        calc2ndderivprobdeltaalpha<-function(x){
          # only to be called inside LLoglmx
          if (x[2]==j & j>1){
            probderiv2deltaalpha<- Z[,x[1]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac0)+frac0*ProbFuncDD(frac0))/vectorsprobs[[j]]
          } else if (x[2]==j+1 & j<no.outcomes){
            probderiv2deltaalpha<-  -Z[,x[1]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac1)+frac1*ProbFuncDD(frac1))/vectorsprobs[[j]]
          } else {
            probderiv2deltaalpha<-vector("numeric",nrow(Z))  
          }
          sum(probderiv2deltaalpha*w)
        }
        
        calc2ndderivprobalphaalpha<-function(x){
          # only to be called inside LLoglmx
          if (x[1]==x[2] & x[1]==j & j>1){
            probderiv2alphaalpha<- -(Zdinv^2)*ProbFuncDD(frac0)/vectorsprobs[[j]]
          } else if (x[1]==x[2] & x[1]==j+1 & j<no.outcomes){
            probderiv2alphaalpha<- (Zdinv^2)*ProbFuncDD(frac1)/vectorsprobs[[j]]
          } else {
            probderiv2alphaalpha<-vector("numeric",nrow(X))
          }
          sum(probderiv2alphaalpha*w)
        }

      hessian<-matrix(vector("numeric",no.parameters^2),nrow=no.parameters,ncol=no.parameters)
      # first add the term that is derived from the cross product of gradients
      if (no.betaparams>0 & no.deltaparams>0){
        scorevecs<-cbind(vectorprobderivbeta,vectorprobderivdelta)
        if (no.threshparams>0){
          scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
        }
        crossprodterms<- -t(scorevecs)%*%(scorevecs*w)
        crossprodterms[upper.tri(crossprodterms)]<-0
        hessian<-hessian+crossprodterms
      } else if (no.betaparams>0){
        scorevecs<-vectorprobderivbeta
        if (no.threshparams>0){
          scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
        }
        crossprodterms<- -t(scorevecs)%*%(scorevecs*w)
        crossprodterms[upper.tri(crossprodterms)]<-0
        hessian<-hessian+crossprodterms
      } else if (no.deltaparams>0){
        scorevecs<-vectorprobderivdelta
        if (no.threshparams>0){
          scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
        }
        crossprodterms<- -t(scorevecs)%*%(scorevecs*w)
        crossprodterms[upper.tri(crossprodterms)]<-0
        hessian<-hessian+crossprodterms
      }
      # then do the same with the cross derivatives
      # clear space in memory, remove first derivatives
      if (no.deltaparams>0){rm(vectorprobderivdelta)}
      if (no.betaparams>0){rm(vectorprobderivbeta)}
      if (no.threshparams>0){rm(vectorprobderivthreshparam)}
      
      if (no.betaparams>0){
        listderivs<-.paircombn(c(1:length(beta))[collectbeta])
        secondderivterms<-sapply(listderivs,calc2ndderivprobbetabeta)
        startrow<-1
        endrow<-no.betaparams
        startcol_1<-0
        counter<-1
        while (counter<=no.betaparams){
          hessian[(startcol_1+counter):endrow,counter+startcol_1]<-hessian[(startcol_1+counter):endrow,counter+startcol_1]+secondderivterms[(no.betaparams*(counter-1)-(counter-2)*(counter-1)/2+1):(no.betaparams*counter-(counter*(counter-1)/2))]
          counter<-counter+1
        }
      }
      
      if (no.deltaparams>0){
        listderivs<-.paircombn(c(1:length(delta))[collectdelta])
        secondderivterms<-sapply(listderivs,calc2ndderivprobdeltadelta)
        startrow<-no.betaparams+1
        endrow<-no.betaparams+no.deltaparams
        startcol_1<-no.betaparams
        counter<-1
        while (counter<=no.deltaparams){
          hessian[(startcol_1+counter):endrow,counter+startcol_1]<-hessian[(startcol_1+counter):endrow,counter+startcol_1]+secondderivterms[(no.deltaparams*(counter-1)-(counter-2)*(counter-1)/2+1):(no.deltaparams*counter-(counter*(counter-1)/2))]
          counter<-counter+1
        }
      }
      
      if (no.threshparams>0){
        listderivs<-.paircombn(c(2:(length(threshparam)-1))[collectthreshparam])
        secondderivterms<-sapply(listderivs,calc2ndderivprobalphaalpha)
        startrow<-no.betaparams+no.deltaparams+1
        endrow<-no.betaparams+no.deltaparams+no.threshparams
        startcol_1<-no.betaparams+no.deltaparams
        counter<-1
        while (counter<=no.threshparams){
          hessian[(startcol_1+counter):endrow,counter+startcol_1]<-hessian[(startcol_1+counter):endrow,counter+startcol_1]+secondderivterms[((no.threshparams*(counter-1)-(counter-2)*(counter-1)/2+1)):(no.threshparams*counter-(counter*(counter-1)/2))]
          counter<-counter+1
        } 
      }
      
      if (no.betaparams>0 & no.deltaparams>0){
        listderivs<-.paircombn(c(1:length(beta))[collectbeta],c(1:length(delta))[collectdelta],same=FALSE)
        secondderivterms<-sapply(listderivs,calc2ndderivprobbetadelta)
        startrow<-no.betaparams+1
        startcol_1<-0
        endrow<-no.betaparams+no.deltaparams
        counter<-1
        while (counter<=no.betaparams){
          hessian[startrow:endrow,counter]<-hessian[startrow:endrow,counter]+secondderivterms[(1+(counter-1)*no.deltaparams):(counter*no.deltaparams)]
          counter<-counter+1
        }
      }
      
      if (no.betaparams>0 & no.threshparams>0){
        listderivs<-.paircombn(c(1:length(beta))[collectbeta],c(2:(length(threshparam)-1))[collectthreshparam],same=FALSE)
        secondderivterms<-sapply(listderivs,calc2ndderivprobbetaalpha)
        startrow<-no.betaparams+no.deltaparams+1
        startcol_1<-0
        endrow<-no.betaparams+no.deltaparams+no.threshparams
        counter<-1
        while (counter<=no.betaparams){
          hessian[startrow:endrow,counter]<-hessian[startrow:endrow,counter]+secondderivterms[(1+(counter-1)*no.threshparams):(counter*no.threshparams)]
          counter<-counter+1
        }
      }
      
      if (no.deltaparams>0 & no.threshparams>0){
        listderivs<-.paircombn(c(1:length(delta))[collectdelta],c(2:(length(threshparam)-1))[collectthreshparam],same=FALSE)
        secondderivterms<-sapply(listderivs,calc2ndderivprobdeltaalpha)
        startrow<-no.betaparams+no.deltaparams+1
        startcol_1<-no.betaparams
        endrow<-no.betaparams+no.deltaparams+no.threshparams
        counter<-1
        while (counter<=no.deltaparams){
          hessian[startrow:endrow,counter+startcol_1]<-hessian[startrow:endrow,counter+startcol_1]+secondderivterms[(1+(counter-1)*no.threshparams):(counter*no.threshparams)]
          counter<-counter+1
        }
      }
  }
  output<-list(scorevector,hessian,BHHHmatrix)  
  }
  if (SameModelMEANSD){
    collectresults<-Map(getLLgradhess,Xr,Xr,weightsr,c(1:no.outcomes))
  } else {
    collectresults<-Map(getLLgradhess,Xr,Zr,weightsr,c(1:no.outcomes))
  }
  scorevector<-Reduce("+",lapply(collectresults,function(x){x[[1]]}))
  
  
  
  
  attr(loglikelihood,"gradient")<-scorevector
  if (analhessian){
    # if coded correctly the hessian calculated up to now is lower triangular, need to make it symmetric
    hessian<-Reduce("+",lapply(collectresults,function(x){x[[2]]}))
    hessian<-hessian+t(hessian)-diag(diag(hessian))
    attr(loglikelihood,"hessian")<-hessian
  }
  if (robustmatrix){
    BHHHmatrix<-Reduce("+",lapply(collectresults, function(x){x[[3]]}))
    attr(loglikelihood,"BHHHhessian")<-BHHHmatrix
  } 
  loglikelihood
}
  # function that calculates the log-likelihood, as a function of the parameters that are not prespecified.
  # the maxLik function used to maximise the log-likelihood requires a function of only the estimated parameters.
  LLoglmxTOP<-function(param){LLoglmx(param,beta=beta,delta=delta,threshparam=threshparam,analhessian=analhessian)}
  # call maxLik to estimate parameters
  maxLikRes<-maxLik(LLoglmxTOP,start=start,iterlim=300,finalHessian=TRUE,method="NR")
  # just remains to extract results and store.
  
  coefficients<-maxLikRes$estimate
  if (robust){
    LLoutput<-LLoglmx(coefficients,beta=beta,delta=delta,threshparam=threshparam,analhessian=analhessian,robustmatrix = TRUE)
    BHHHmatrix<-attr(LLoutput,"BHHHhessian")
  } else {
    BHHHmatrix<-NULL
  }
  
  # separate coefficients into beta, delta and alpha, with prespecified values
  betacoeffs<-deltacoeffs<-threshparamcoeffs<-vector("logical",length(coefficients))
  if (no.betaparams>0){
    beta[is.na(beta)]<-coefficients[1:no.betaparams]
    betacoeffs[1:no.betaparams]<-TRUE
  }
  if (no.deltaparams>0){
    delta[is.na(delta)]<-coefficients[(no.betaparams+1):(no.betaparams+no.deltaparams)]
    deltacoeffs[(no.betaparams+1):(no.betaparams+no.deltaparams)]<-TRUE
  }
  if (no.threshparams>0){
    threshparam[is.na(threshparam)]<-coefficients[(no.betaparams+no.deltaparams+1):(no.betaparams+no.deltaparams+no.threshparams)]
    threshparamcoeffs[(no.betaparams+no.deltaparams+1):(no.betaparams+no.deltaparams+no.threshparams)]<-TRUE
  }
  allparameters<-list(beta,delta,threshparam)
  
  names(allparameters)<-c("beta","delta","threshparam")
  coeff_type<-list(betacoeffs,deltacoeffs,threshparamcoeffs)
  names(coefficients)<-outputnames
  attr(coefficients,"coefftypes")<-coeff_type
  loglikelihood<-maxLikRes$maximum
  #attr(loglikelihood,"BaselineLL")<-BaselineLL
  attr(loglikelihood,"No.Obs")<-No.Obs
  if (weighted){
    weights<-weights
  } else {
    weights<-NULL
  }

  results<-list(loglikelihood=loglikelihood,link=link,no.iterations=maxLikRes$iterations,coefficients=coefficients,returnCode=maxLikRes$code,call=call,gradient=maxLikRes$gradient,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData
                ,hessian=maxLikRes$hessian,BHHHhessian=BHHHmatrix,Hetero=Heteroskedastic,NOutcomes=no.outcomes,Outcomes=listoutcomes,BothEq=BothMeanVar,sdmodel=sdmodel,allparams=allparameters,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary),Est.Parameters=Est.Parameters,modelframes=modelframes)
  class(results)<-c("oglmx")
  invisible(results)
}

vcov.oglmx<-function(object,tol=1e-20,...){
  if (is.null(object$BHHHhessian)){
    vcov<-qr.solve(-object$hessian,tol=tol)
  } else {
    vcov<- qr.solve(object$hessian,tol=tol)%*%(object$BHHHhessian*(attr(object$loglikelihood,"No.Obs")/(attr(object$loglikelihood,"No.Obs")-1)))%*%qr.solve(object$hessian,tol=tol)
  }
  colnames(vcov)<-rownames(vcov)<-names(object$coefficients)
  return(vcov)
}

summary.oglmx<-function(object,tol=1e-20, ... ){
  stdEr.oglmx<-diag(vcov(object,tol=tol))^0.5
  t<-object$coefficients/stdEr.oglmx
  p <- 2*pnorm( -abs( t))
  results <- cbind("Estimate"=object$coefficients,
                   "Std. error"=stdEr.oglmx,
                   "t value"=t, "Pr(>|t|)"=p)
  betaresults<-results[attr(object$coefficients,"coefftypes")[[1]], ,drop=FALSE]
  deltaresults<-results[attr(object$coefficients,"coefftypes")[[2]], ,drop=FALSE]
  cutoffresults<-results[attr(object$coefficients,"coefftypes")[[3]], ,drop=FALSE]
  resultsSplit<-list(betaresults,deltaresults,cutoffresults)
  summary<-list(regtype=.regtype.oglmx(object),loglikelihood=object$loglikelihood,estimate=results,estimateDisplay=resultsSplit,no.iterations=object$no.iterations,McFaddensR2=McFaddensR2.oglmx(object),AIC=AIC(object),coefficients=object$coefficients)
  class(summary)<-"summary.oglmx"
  summary
}

print.summary.oglmx<-function(x, ... ){
  cat(x$regtype,"\n")
  cat("Log-Likelihood:", x$loglikelihood, "\n")
  cat("No. Iterations:", x$no.iterations, "\n")
  cat("McFadden's R2:",x$McFaddensR2,"\n")
  cat("AIC:",x$AIC,"\n")
  if (nrow(x$estimateDisplay[[1]])>0 & nrow(x$estimateDisplay[[2]])==0 & nrow(x$estimateDisplay[[3]])==0){
    printCoefmat(x$estimateDisplay[[1]])
  } else if (nrow(x$estimateDisplay[[1]])>0){
    if (nrow(x$estimateDisplay[[2]])>0){
      cat("-----","Mean Equation","------\n")
    }
    printCoefmat(x$estimateDisplay[[1]],signif.legend=FALSE)
  }
  if (nrow(x$estimateDisplay[[2]])>0){
    if (nrow(x$estimateDisplay[[1]])>0){
      cat("-----","SD Equation","------\n")
    }
    if (nrow(x$estimateDisplay[[3]])>0){
      printCoefmat(x$estimateDisplay[[2]],signif.legend=FALSE)
    } else {
      printCoefmat(x$estimateDisplay[[2]])
    }
  }
  if (nrow(x$estimateDisplay[[3]])>0){
    cat("-----","Threshold Parameters","-----\n")
    printCoefmat(x$estimateDisplay[[3]])
  }
}

.BaseLL<-function(object){
  #outcome<-peersimdata$score
  #weight<-(1/peersimdata$rankprob)*5760/sum((1/peersimdata$rankprob))
  #values<-as.numeric(levels(as.factor(outcome)))
  #countoutcomes<-vector("numeric",0)
  #for (i in 1:length(values)){
  #  countoutcomes[i]<-sum(as.numeric(outcome)==values[i])
  #}
  #if (sum(weight == rep_len(1,length(outcome)))==length(outcome)){
  #  BaseLL<-sum(countoutcomes*log(countoutcomes/length(outcome)))
  #} else {
  #  BaseLL<-0
  #  for (i in 1:length(values)){
  #    BaseLL<-BaseLL+sum(weight[as.numeric(outcome)==values[i]]*log(countoutcomes[i]/length(outcome)))
  #  }
  #}
 BaseLL<-as.numeric(logLik(oglmx(object$NoVarModData$Y~1,data=object$NoVarModData,weights=object$NoVarModData$weights)))
 return(BaseLL)
}

McFaddensR2.oglmx<-function(object){
  value<-1-logLik(object)/.BaseLL(object)
  return(value)
}

AIC.oglmx<-function(object, ..., k=2){
  # 2*number of estimatated parameters - 2*log likelihood
  value<-k*length(object$coefficients)-2*logLik(object)
  return(value)
}

logLik.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

logLik.summary.oglmx<-function(object, ...){
  object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

coef.oglmx<-function(object, ...){
  coefnames<-names(object$coefficients)
  output<-as.vector(object$coefficients)
  names(output)<-coefnames
  return(output)
}

coef.summary.oglmx<-function(object, ...){
  coefnames<-names(object$coefficients)
  output<-as.vector(object$coefficients)
  names(output)<-coefnames
  return(output) 
}

nobs.oglmx<-function(object, ...){
  return(attr(object$loglikelihood,"No.Obs"))
}

formula.oglmx<-function(x, ...){
  # extract the formula for an oglmx object
  # for use to apply a model name in lrtest
  if (is.null(x$formula[[2]])){
    value<-x$formula[[1]]
  } else {
    # collect the names from the terms output
    # from the mean equation, include response term
    meannames<-names(attr(terms(x)[[1]],"dataClasses"))
    varnames<-attr(terms(x)[[2]],"term.labels")
    textoutput<-paste(meannames[1],"~",meannames[2])
    if (length(meannames)>2){
      for (j in 3:length(meannames)){
        textoutput<-paste(textoutput,"+",meannames[j])
      }
    }
    textoutput<-paste(textoutput,"|",varnames[1])
    if (length(varnames)>1){
      for (j in 2:length(varnames)){
        textoutput<-paste(textoutput,"+",varnames[j])
      }
    }
    value<-formula(textoutput)
  }
  return(value)
}

.regtype.oglmx<-function(object){
  if (sum(object$NoVarModData$weights==1)!=nrow(object$NoVarModData)){
    Zero<-"Weighted "
  } else {
    Zero<-""
  }
  if (object$Hetero){
    First<-"Heteroskedastic "
  } else {
    First<-""
  }  
  if (object$NOutcomes>2){
    Second<-"Ordered "
  } else {
    Second<-""
  } 
  if (object$link=="logit"){
    Third<-"Logit "
  } else if (object$link=="probit"){
    Third<-"Probit "
  } else if (object$link=="cloglog"){
    Third<-"CLogLog "
  } else if (object$link=="loglog"){
    Third<-"LogLog "
  } else if (object$link=="cauchit"){
    Third<-"Cauchit "
  }
  Fourth<-"Regression"
  value<-paste(Zero,First,Second,Third,Fourth,sep="")
  return(value)
}

.paircombn<-function(vec1,vec2=NULL,same=TRUE){
  result1<-vector("numeric",0)
  result2<-vector("numeric",0)
  if (same){
    for (i in 1:length(vec1)){
      result1<-c(result1,rep(vec1[i],length(vec1)+1-i))
      result2<-c(result2,vec1[i:length(vec1)])           
    }
  } else {
    if (is.null(vec2)){stop("If pairs are not drawn from the same list argument vec2 should be specified.")}
    for (i in 1:length(vec1)){
      result1<-c(result1,rep(vec1[i],length(vec2)))
      result2<-c(result2,vec2)
    }
  }
  lapply(1:length(result1),function(x){c(result1[x],result2[x])})
}

.checkbinary<-function(x){
  if (sum(x==0)+sum(x==1)==length(x) & sum(x==0)>0 & sum(x==1)>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# function to check if there is an NA in the row
.IsNARow<-function(x){
  if (sum(is.na(x))>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

probit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,analhessian=TRUE,na.action=TRUE,savemodelframe=FALSE,robust=FALSE){
  value<-oglmx(formulaMEAN=formula,data=data,start=start,beta=beta,analhessian=analhessian,na.action=na.action,savemodelframe=savemodelframe,link="probit",constantMEAN=TRUE,constantSD=FALSE,delta=0,threshparam=0)
  return(value)
}

oprobit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,analhessian=TRUE,na.action=TRUE,savemodelframe=FALSE,robust=FALSE){
  value<-oglmx(formulaMEAN=formula,data=data,start=start,beta=beta,analhessian=analhessian,na.action=na.action,savemodelframe=savemodelframe,link="probit",constantMEAN=FALSE,constantSD=FALSE,delta=0,threshparam=NULL,robust=robust)
  return(value)
}

ologit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,analhessian=TRUE,na.action=TRUE,savemodelframe=FALSE,robust=FALSE){
  value<-oglmx(formulaMEAN=formula,data=data,start=start,beta=beta,analhessian=analhessian,na.action=na.action,savemodelframe=savemodelframe,link="logit",constantMEAN=FALSE,constantSD=FALSE,delta=0,threshparam=NULL,robust=FALSE)
  return(value)
}

logit.reg<-function(formula,data,start=NULL,weights=NULL,beta=NULL,analhessian=TRUE,na.action=TRUE,savemodelframe=FALSE,robust=FALSE){
  value<-oglmx(formulaMEAN=formula,data=data,start=start,beta=beta,analhessian=analhessian,na.action=na.action,savemodelframe=savemodelframe,link="logit",constantMEAN=TRUE,constantSD=FALSE,delta=0,threshparam=0,robust=FALSE)
  return(value)
}

