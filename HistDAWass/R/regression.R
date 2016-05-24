## Regression analysis ----
# Multiple regression analysis of histogram variable based on Wasserstein distance----
#' Multiple regression analysis for  histogram variables based on a two component model and L2 Wasserstein distance
#' @description The function implements Multiple regression analysis for  histogram variables based on 
#' a two component model and L2 Wasserstein distance. Taking as imput dependent histogram variable and
#'  a set of explanatory histogram variables the methods return a least squares estimation of a two component
#'  regression model based on the decomposition of L2 Wasserstein metric for distributional data.
#' 
#' @param data A MatH object (a matrix of distributionH).
#' @param Yvar An integer, the dependent variable number in data.
#' @param Xvars A set of integers the explanantory variables in data.
#' @param simplify a logical argument (default=FALSE). If TRUE only few equally spaced quantiles 
#' are considered (for speeding up the  algorithm)
#' @param qua If \code{simplify=TRUE} is the number of quantiles to consider.

#' @return a named vector with the model estimated parameters
#' @references 
#' Irpino A, Verde R (in press 2015). Linear regression for numeric symbolic variables: a least squares approach 
#' based on Wasserstein Distance. ADVANCES IN DATA ANALYSIS AND CLASSIFICATION, ISSN: 1862-5347, DOI:10.1007/s11634-015-0197-7 \cr
#' An extended version is available  on arXiv repository arXiv:1202.1436v2 \url{http://arxiv.org/abs/1202.1436v2}
#' @details  A two component regression model is implemented. The observed variables are histogram variables
#'  according to the definition given in the framework 
#' of Symbolic Data Analysis and the parameters of the model are estimated 
#' using the classic Least Squares method. An appropriate metric is introduced 
#' in order to measure the error between the observed and the predicted distributions. 
#' In particular, the Wasserstein distance is proposed. 
#' Such a metric permits to predict the response variable as direct linear combination of other independent
#'  histogram variables.
#' @examples
#' model.parameters=WH.regression.two.components(data = BLOOD,Yvar = 1, Xvars= c(2:3))
#' @importFrom stats as.formula lm
#' @export
WH.regression.two.components=function(data,Yvar,Xvars,simplify=FALSE,qua=20){
  #check input
  if (is(data)!="MatH"){stop("HW.regression.two.components accepts only MatH objects as input data table")}
  if(length(Yvar)>1){stop("Multiple choice not allowed for Y variable")}
  #check for missing values and set up working data
  selected=c(1:nrow(data@M))
  Y=data[selected,Yvar]
  X=data[selected,Xvars]
  n=length(selected)
  d=ncol(X@M)
  #simplify, i.e. use a fixed number of quantiles for a rapid computation
  
  if (simplify){
    pr=c(0:qua)/qua
    for (i in 1:n){
      tmpx=numeric(0)
      for (q in 0:qua){
        tmpx=c(tmpx, compQ(Y@M[i,1][[1]],pr[q+1]))
      }
      Y@M[i,1][[1]]=new("distributionH",x=tmpx,p=pr)
      for (j in 1:d){
        tmpx=numeric(0)
        for (q in 0:qua){
          tmpx=c(tmpx, compQ(X@M[i,j][[1]],pr[q+1]))
        }
        X@M[i,j][[1]]=new("distributionH",x=tmpx,p=pr)
      }
    }
  }
  #
  #extract means and do multiple regression on means
  MatAver=matrix(0,n,(d+1))#the matrix of means of the distributions
  for (i in 1:n){
    MatAver[i,1]=Y@M[i,1][[1]]@m
    for (j in 1:d){
      MatAver[i,(j+1)]=X@M[i,j][[1]]@m
    }
  }
  colnames(MatAver)=paste0("AV_",c(colnames(Y@M), colnames(X@M)))
  rownames(MatAver)=rownames(Y@M)
  MatAver=as.data.frame(MatAver)
  xnam = colnames(MatAver)[2:(d+1)]
  fmla = as.formula(paste(colnames(MatAver)[1], paste(xnam, collapse= "+"),sep="~"))
  fit = lm(fmla, data=MatAver)
  AveCoeff=fit$coefficient
  names(AveCoeff)[1]="(AV_Intercept)"
  
  #center data
  CENDATA=new("MatH",nrows=n,ncols=d+1)
  for (i in 1:n){
    CENDATA@M[i,1][[1]]=Y@M[i,1][[1]]-Y@M[i,1][[1]]@m
    for(j in 1:d){
      CENDATA@M[i,(j+1)][[1]]=X@M[i,j][[1]]-X@M[i,j][[1]]@m  
    }
  }
  if (!simplify){
    CENDATA=registerMH(CENDATA)
    
    #register all data
  }
  YC=CENDATA[,1]
  XC=CENDATA[,2:(d+1)]
  #do Non Negative Least Squares modifying
  # Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  gammas=WH.NNLS(XC,YC)
  gammas=as.vector(gammas)
  names(gammas)=paste0("CEN_",colnames(X@M))
  return(parameters=c(AveCoeff,gammas))
}
## NNLS for histogram variables -----

WH.NNLS=function(X,Y){
  # Non Negative Least Squares for histogram variables modifying
  # Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  tol=1e-14
  #Check input
  if ((is(Y)[1]!=is(X)[1])&&(is(X)[1]!="MatH")){stop("X and Y maust be a MatH objects")}  
  if(nrow(Y@M)!=nrow(X@M)){stop("X and Y must have the same number of rows")}
  #Lawson and Hanson NNLS algorithm
  # step 1
  P=integer(0)
  Z=c(1:ncol(X@M))
  gamma=matrix(0,nrow = ncol(X@M),ncol = 1)
  stp=0
  # step 2
  while (stp==0){
    
    w=WH.mat.prod(X,Y,traspose1=TRUE)-WH.mat.prod(X,X,traspose1=TRUE)%*%gamma
    
    #step 3: if Z is empty or if all w are less or equal to 0 end
    if ((length(Z)==0)||(sum(w<=0)==length(w))){
      stp=1
    }
    else{
      P=c(P,Z[which.max(w[Z,1])])
      Z=Z[-which.max(w[Z,1])]
      stp2=0
      while (stp2==0){
        P=sort(P)
        Z=sort(Z)
        tmpX=X[,P]
        GG=solve(WH.mat.prod(tmpX,tmpX,traspose1=TRUE))%*%WH.mat.prod(tmpX,Y,traspose1=TRUE)
        zvect=matrix(0,ncol(X@M),1)
        zvect[P,1]=GG
        #step 7
        if (length(zvect[P,1])==sum(zvect[P,1]>0)){
          gamma=zvect
          stp=0
          stp2=1
          #go to step 2
        }
        else{#step 8 and 9
          tmp3=gamma/(gamma-zvect)
          alpha=min(tmp3[zvect[P,1]<=0,1])
          #step 10
          gamma=gamma+alpha*(zvect-gamma)
          
          #step 11
          #remove from P all the indices for wich gamma_j is 0 and move them to Z
          tmp4=P[abs(gamma[P,1])<tol]
          P=P[-tmp4]
          Z=c(Z,tmp4)
          }
      }
    }
  }
  return(gamma) 
}
## diagnostics measures to be implemented
#' Multiple regression analysis for  histogram variables based on a two component model and L2 Wasserstein distance
#' @description Predict distributions using the results of a regression done with \code{WH.regression.two.components} function.
#' 
#' @param data A MatH object (a matrix of distributionH) explantory part.
#' @param parameters A named vector with the parameter from a \code{WH.regression.two.components} model
#' @return a \code{MatH}  object, the predicted histograms
#' @references 
#' Irpino A, Verde R (in press 2015). Linear regression for numeric symbolic variables: a least squares approach 
#' based on Wasserstein Distance. ADVANCES IN DATA ANALYSIS AND CLASSIFICATION, ISSN: 1862-5347, DOI:10.1007/s11634-015-0197-7 \cr
#' An extended version is available  on arXiv repository arXiv:1202.1436v2 \url{http://arxiv.org/abs/1202.1436v2}
#' @examples
#' # do regression
#'  model.parameters=WH.regression.two.components(data = BLOOD,Yvar = 1, Xvars= c(2:3))
#' # do prediction
#' Predicted.BLOOD=WH.regression.two.components.predict(data = BLOOD[,2:3],parameters=model.parameters)
#' @export
WH.regression.two.components.predict=function(data,parameters){
  if (is(data)!="MatH"){stop("HW.regression.two.components.predict accepts only MatH objects as input data table")}
  if (ncol(data@M)!=(length(parameters)-1)/2){stop("Predictors (colums of data) are not consistent with the number of parameters (that must be (num.of colums of data)x2+1)")}
  npred=ncol(data@M);
  if (sum(abs(parameters[(1+npred+1):length(parameters)])-(parameters[(1+npred+1):length(parameters)]))>0){
    stop(cat("HW.regression.two.components.predict accepts only last ",npred," parameters positives"))
  }
  indiv=nrow(data@M)
  predictions=new("MatH",nrows=indiv, ncols=1, names.cols=c("Predicted"), names.rows=rownames(data@M))
  parameters=as.numeric(parameters)
  for (i in 1:indiv){
    tmp.pred=new("distributionH",x=c(parameters[1],parameters[1]),p=c(0,1))
    for (j in 1:npred){
      tmpdist=data@M[i,j][[1]]
      tmp.pred=tmp.pred+(parameters[1+j]*tmpdist@m)+(parameters[1+npred+j]*(tmpdist-tmpdist@m))
      
    }
    predictions@M[i,1][[1]]=tmp.pred
  }
  return(predictions)
}
#' Goodness of Fit indices for Multiple regression of histogram variables based on a two component model and L2 Wasserstein distance
#' @description It computes three goodness of fit indices using the results and the predictions of a regression done with \code{WH.regression.two.components} function.
#' @param observed A one column MatH object, the observed histogram variable
#' @param predicted A one column MatH object, the predicted histogram variable.
#' @return a list with the GOF indices
#' @references 
#' Irpino A, Verde R (in press 2015). Linear regression for numeric symbolic variables: a least squares approach 
#' based on Wasserstein Distance. ADVANCES IN DATA ANALYSIS AND CLASSIFICATION, ISSN: 1862-5347, DOI:10.1007/s11634-015-0197-7 \cr
#' An extended version is available  on arXiv repository arXiv:1202.1436v2 \url{http://arxiv.org/abs/1202.1436v2}
#' @examples
#' # do regression
#'  model.parameters=WH.regression.two.components(data = BLOOD,Yvar = 1, Xvars= c(2:3))
#'  #' # do prediction
#' Predicted.BLOOD=WH.regression.two.components.predict(data = BLOOD[,2:3],parameters=model.parameters)
#' # compute GOF indices
#' GOF.indices=WH.regression.GOF(observed=BLOOD[,1], predicted=Predicted.BLOOD)
#' @export
WH.regression.GOF=function(observed, predicted){
    indiv=nrow(observed@M)
  WassD2=matrix(0,nrow = indiv,ncol = 1)
  TOTSSQ=WH.SSQ(observed)
  MO=WH.vec.mean(observed)
  SSQR=0
  RMSE_W=0;
  NUM_OMEGA=0
  DEN_OMEGA=0
  for (i in 1:indiv){
    WassD2[i,1]=WassSqDistH(observed@M[i,1][[1]],predicted@M[i,1][[1]])
    SSQR=SSQR+WassSqDistH(predicted@M[i,1][[1]], MO)
    NUM_OMEGA=NUM_OMEGA+(predicted@M[i,1][[1]]@m-MO@m)^2+(predicted@M[i,1][[1]]@s)^2
    DEN_OMEGA=DEN_OMEGA+(observed@M[i,1][[1]]@m-MO@m)^2+(observed@M[i,1][[1]]@s)^2
  }

  return(indices=list(RMSE_W=sqrt(sum(WassD2)/indiv), 
                      OMEGA=NUM_OMEGA/DEN_OMEGA,
                      PSEUDOR2=list(index=max(0, min(SSQR/TOTSSQ,1-sum(WassD2)/TOTSSQ)),
                                    details=c(TotSSQ=TOTSSQ, SSQ.R=SSQR, SSQ.E=sum(WassD2), 
                                              Bias=TOTSSQ-SSQR-sum(WassD2), SSQ.R.rel=SSQR/TOTSSQ,
                                              SSQ.E.rel=sum(WassD2)/TOTSSQ, 
                                              SSQ.bias.rel=(TOTSSQ-SSQR-sum(WassD2))/TOTSSQ))))
}

WH.regr.two.comp.Bootstrap=function(data,Yvar,Xvars,simplify=FALSE,
                                    qua=20,rep=100, GOF=FALSE){
  nobj=get.MatH.nrows(data)
  ind=sample(nobj,nobj,replace = TRUE)
  pars=WH.regression.two.components(data[ind,],Yvar,Xvars,simplify=simplify,qua=qua)
  
  if (GOF==TRUE){
    OBS=data[ind,Yvar]
    PRED=WH.regression.two.components.predict(data[ind,Xvars],pars)
    idxs=WH.regression.GOF(OBS,PRED)
    MP=c(pars,RMSEW=idxs$RMSE_W, OMEGA=idxs$OMEGA,PSEUDOR2=idxs$PSEUDOR2$index)
  }else{
    MP=pars
  }
  for (i in 1:rep){
    ind=sample(nobj,nobj,replace = TRUE)
    pars=WH.regression.two.components(data[ind,],Yvar,Xvars,simplify=simplify,qua=qua)
    
    if (GOF==TRUE){
    OBS=data[ind,Yvar]
    PRED=WH.regression.two.components.predict(data[ind,Xvars],pars)
    idxs=WH.regression.GOF(OBS,PRED)
    TMP=c(pars,RMSEW=idxs$RMSE_W, OMEGA=idxs$OMEGA,PSEUDOR2=idxs$PSEUDOR2$index)
    MP=rbind(MP,TMP)
    }else{
      MP=rbind(MP,pars)
    }
  }
  print(summary(MP))
  return(MP)
}
