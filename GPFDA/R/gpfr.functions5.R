main1=function(response,lReg=NULL,fReg=NULL,fxList=NULL,fbetaList=NULL){
  ### main1 is the function for the case that the response variable is 
  ### one dimensional. This function will do lm, flm, gp, or conbinations 
  ### of the three, depends on the avaliable data.
  y=response
  
  ## multivariate linear regression
  ml=NULL
  if (!is.null(lReg)){
     #  if(class(lReg)!='matrix') stop('The covariates for scaler multivariate linear regression are expected to store in a matrix, not other format')
    ml=lm(y~lReg)
    resid_ml=as.matrix(y-resid(ml))
    y=resid_ml
  }
  
  ## functional regression
  temp=list(NULL)
  if(!is.null(fReg)){
    if(class(fReg)=='matrix'|class(fReg)=='fd')
      fReg=list(fReg)
    if(class(fReg)=='list'){
      if(length(unique(unlist(lapply(fReg,class))))!=1) 
        stop('functional covariates are expected to have same class')
      if(unique(unlist(lapply(fReg,class)))=='matrix'){
        temp=list(NULL)
        res=list(y)
        
        ## set up functional variable for fx
        if(length(fxList)!=length(fReg)){
          cat('Length of fx list is not equal to the length of list of functional covariates','\n')
          
          if(length(fxList)==0){
            cat('     Defualt fx list is applied','\n')
            fxList=lapply(fReg,function(i){
              i=min(as.integer(ncol(i)/5),10)
              names(i)=list('nx_basis')
              return(i)
            })
          } 
          else{
            cat('     First item in fx list is applied to all items','\n')
            fxList=lapply(fReg,function(i){
              i=fxList[[1]]
              return(i)
            })
          }
          
        }
          
        ## set up functional parameter for fbeta
        if(length(fbetaList)!=length(fReg)){
          cat('Length of fbeta list is not equal to the length of list of functional covariates','\n')
          if(length(fbetaList)==0){
            cat('     Defualt fbeta list is applied','\n')
            fbetaList=fxList
            fbetaList=lapply(fbetaList,function(i){
              names(i)='nbeta_basis'
              return(i)
            })
          }
          else{
            cat('     First item in fbeta list is applied to all items','\n')
            fbetaList=lapply(fReg,function(i){
              i=fbetaList[[1]]
              return(i)
            })
          }
        }
          
        for(i in seq_along(fReg)){
          ## functional regression with scaler response and functional covariates using fdata.
          mf=freg1(y,fReg[[i]],fxList[[i]],fbetaList[[i]])
          res=list(res,as.matrix(resid(mf)))
          temp=c(temp,list(mf))
          if(is.null(temp[[1]])) temp=temp[-1]
          y=res[[length(res)]]
        } 
      }
      if(unique(unlist(lapply(fReg,class)))=='fd'){
        res=list(y)
        ## set up functional parameter for fbeta
        if(length(fbetaList)!=length(fReg)){
          cat('     Length of fbeta list is not equal to the length of list of functional covariates','\n')
          if(length(fbetaList)==0){
            cat('     Defualt fbeta list is applied','\n')
            fbeta[[i]]=betaPar()
          }
          else{
            cat('     First item in fbeta list is applied to all items','\n')
            fbeta=lapply(fReg,function(i){
              i=betaPar(fbetaList[[1]])
              return(i)
            })
          }
        }
        if(length(fbetaList)==length(fReg))
          fbeta=lapply(fbetaList,betaPar)
          
        for(i in seq_along(fReg)){          
          ## functional regression with scaler response and functional covariates using fd.
          mf=freg2(y,fReg[[i]],fbeta[[i]])
          res=list(res,as.matrix(resid(mf)))
          temp=c(temp,list(mf))
          if(is.null(temp[[1]])) temp=temp[-1]
          y=res[[length(res)]]
        } 
      }
    }
  }
  out=list(model=list(ml=ml,mf=temp),res=y)
  return(out)
}

freg1=function(y,fx,nx_bas,nbeta_bas){
  ### works for scaler response and single functional covariates, 
  ### return the regression model.
  ### y is expected to be a column matrix; fx is expected to be a matrix
  if (nrow(fx)!=nrow(y)) stop('unequal sample size for response and functional covariates')
  fx=fdata(fx,argvals=seq(0,1,len=ncol(fx)))
  tt=fx[['argvals']]
  basis1=create.bspline.basis(rangeval=range(tt),nbasis=nx_bas)
  basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbeta_bas)
  fm=fregre.basis(fx,y,basis1,basis2)
  return(fm)
}

freg2=function(y,xlist,fbeta){
  ### works for scaler response and single functional covariates,
  ### return the regression model
  ### y is expected to be a column matrix; fx is expected to be an fd object
  ### fbeta is expected to be an fd onject
  fm=fRegress(y[,1], xlist, fbeta)
  return(fm)
}



main2=function(response,lReg=NULL,fReg=NULL,fyList=NULL,fbetaList_l=NULL,
               fxList=NULL,concurrent=TRUE,fbetaList_f=NULL,time=NULL){
  ### main2 is for the regression with functional response, 
  ### This function will do lm, flm, gp, or conbinations 
  ### of the three, depends on the avaliable data.
  y=response
  ml=NULL;res=NULL;fittedFM=NULL;
  if(!is.null(lReg)){
    if(class(y)=='fdata') y=y$data
  }
  if(!is.null(time)){
    if(!is.null(fyList)) fyList$time=time
    if(is.null(fyList)) fyList=list(time=time)
    if(!is.null(fbetaList_l)) fbetaList_l=lapply(fbetaList_l,function(i) c(i,list(rtime=range(time))))
    if(is.null(fbetaList_l)) fbetaList_l=list(list(rtime=range(time)))
    if(!is.null(fbetaList_f)) fbetaList_f=lapply(fbetaList_f,function(i) c(i,list(rtime=range(time))))
    if(is.null(fbetaList_f)) fbetaList_f=list(list(rtime=range(time)))
    if(!is.null(fxList)) fxList=lapply(fxList,function(i) c(i,list(rtime=range(time))))
    if(is.null(fxList)) fxList=list(list(time=time))
  }
  
  if(class(y)=='matrix'){
    ## define fd object for y if y is a matrix
    y=mat2fd(y,fyList)
  }
  if(class(y)!='fd') stop('class of response must be one of matrix, fd or fdata')
  y_time=seq(y$basis$rangeval[1],y$basis$rangeval[2],len=length(y$fdnames$time))
  
  if(is.null(fbetaList_l[[1]]$nbasis)) fbetaList_l=lapply(fbetaList_l,function(i)c(i,list(nbasis=y$basis$nbasis)))
  if(is.null(fbetaList_l[[1]]$norder)) fbetaList_l=lapply(fbetaList_l,function(i)c(i,list(norder=c(fyList$norder,6)[1])))
  if(is.null(fbetaList_l[[1]]$Pen)) {fbetaList_l=lapply(fbetaList_l,function(i){
    if(!is.null(fyList$Pen)) c(i,list(Pen=fyList$Pen))
    if(is.null(fyList$Pen)) c(i,Pen=c(0,0))
  })}
    
  if(!is.null(lReg)){
    ## define list of x 
    if(class(lReg)!='matrix') stop('class of lReg is expected to be matrix')
    x=lReg
    nx=ncol(x)
    lxList=vector('list',length=nx)
    for(i in 1:nx) lxList[[i]]=x[,i]
    
    ## define list of beta
    if(length(fbetaList_l)!=length(lxList)){
      cat('     Length of fbetaList_l list is not equal to the length of list of functional covariates','\n')
      if(length(fbetaList_l)==0){
        cat('     Default fbetaList_l is applied ','\n')
        betalist=lapply(lxList,function(i){
          i=betaPar()
        })
      }
      
      if(length(fbetaList_l)>0){
        cat('    The first fbetaList_l is applied to all items', '\n')
        betalist=lapply(lxList,function(i){
          i=betaPar(fbetaList_l[[1]])
        })
      }
    }
    if(length(fbetaList_l)==length(lxList))
      betalist=lapply(fbetaList_l,betaPar)
    
    
    #regression
    ml = fRegress(y, lxList, betalist)
    betaEstMat=do.call('cbind',lapply(ml$betaestlist,function(i) predict(i,y_time)))
    ml_fitted=lReg%*%t(betaEstMat)
    
    if(class(response)=='fd') y_raw=eval.fd(y_time,response)
    if(class(response)=='matrix'){
      if(nrow(response)==nrow(ml_fitted)) residML=response-ml_fitted
      if(nrow(response)==ncol(ml_fitted)) residML=t(response)-ml_fitted
    }
    y=residML
    res=c(res,list(y));fittedFM=c(fittedFM,list(ml_fitted))
  }
  mfTrainfd=NULL
  ## functional response with functional covariates
  temp=list(NULL)
  if(!is.null(fReg)){
    y=mat2fd(y,fyList)
    
    ## set up list of fd object for x
    if(class(fReg)=='matrix' | class(fReg)=='fd')
      fReg=list(fReg)
    if(class(fReg)=='list'){
      if(length(unique(unlist(lapply(fReg,class))))!=1) 
        stop('functional covariates are expected to have same class')
      if(unique(unlist(lapply(fReg,class)))=='matrix'){
        if(ncol(fReg[[1]])!=length(y$fdnames$time)) fReg=lapply(fReg,t)
        if(length(fxList)!=length(fReg)){
          cat('     Length of fxList list is not equal to the length of list of functional covariates','\n')
          fReg=lapply(fReg,t)
          if(length(fxList)==0){
            cat('     Default fxList is applied','\n')
            fReg=lapply(fReg,mat2fd)
          }
          if(length(fxList)>0){
            cat('     First fxList is applied to all items','\n')
            
            fReg=lapply(fReg,function(i){
              i=mat2fd(i,fxList[[1]])
            })
          }
        }
        if(length(fxList)==length(fReg)){
          fReg=lapply(1:length(fReg),function(i){
            mat2fd(fReg[[i]],fxList[[i]])
          })
        }
      }
      
      
      ## set up list of fdPar object for beta
      if(length(fbetaList_f)!=length(fReg)){
        cat('     Length of fbetaList_f list is not equal to the length of list of functional covariates','\n')
        if(length(fbetaList_f)==0){
          cat('     Default fbetaList_f is applied','\n')
          if(1-concurrent) betalist=lapply(fReg,function(i) i=betaPar(list(bivar=TRUE)))
          if(concurrent) betalist=lapply(fReg,function(i) i=betaPar())
        }
        if(length(fbetaList_f)>0){
          cat('     First fbetaList_f is applied to all items','\n')
          if(1-concurrent) betalist=lapply(fReg,function(i) i=betaPar(list(bivar=TRUE)))
          if(concurrent) betalist=lapply(fReg,function(i) i=betaPar(fbetaList_f[[1]]))
        }
      }
      if(length(fbetaList_f)==length(fReg)){
        if(1-concurrent) betalist=lapply(fbetaList_f, function(i) betaPar(list(bivar=TRUE)))
        if(concurrent) betalist=lapply(fbetaList_f, betaPar)
      }
        
      ## regression
      temp=NULL
      for(i in seq_along(fReg)){
        if(is.matrix(y)) y=mat2fd(y,fyList)
        x=fReg[[i]]
        if(concurrent){
          const=rep(1,dim(x$coef)[2])
          xlist=list(const=const,x=x)
          
          b1=betalist[[i]]
          bList=list(const=b1,x=b1)
          mf=fRegress(y,xlist,bList) 
          temp=c(temp,list(mf))
          
          # evaluate functional coefficients
          betaEstMat=list(predict(const=mf$betaestlist$const,y_time)[,1],x=predict(const=mf$betaestlist$x,y_time)[,1])
          mf_fitted=apply(t(eval.fd(y_time,x)),1,function(i) i=i*betaEstMat[[2]]+betaEstMat[[1]])
          
          if(class(y)=='fd') y_raw=t(eval.fd(y_time,y))
          if(class(y)=='matrix'){
            if(nrow(y)==nrow(mf_fitted)) residMF=y_raw-mf_fitted
            if(nrow(y)==ncol(mf_fitted)) residMF=t(y_raw)-mf_fitted
          }
          y=residMF
          res=c(res,list(y));fittedFM=c(fittedFM,list(mf_fitted))
        }
        
        if(1-concurrent){
          bList = list(betaPar(), betalist[[i]])
          mf=linmod(x,y,bList) 
          temp=c(temp,list(mf))
          
          betaEstMat=list(b0=eval.fd(y_time,mf$beta0estfd)[,1],b1=eval.bifd(y_time,y_time,mf$beta1estbifd)[,1])
          mf_fitted=apply(t(eval.fd(y_time,x)),1,function(i) i=i%*%betaEstMat[[2]]/length(y_time)^2+betaEstMat[[1]])
          
          if(class(y)=='fd') y_raw=t(eval.fd(y_time,y))
          if(class(y)=='matrix'){
            if(nrow(y)==nrow(mf_fitted)) residMF=y_raw-mf_fitted
            if(nrow(y)==ncol(mf_fitted)) residMF=t(y_raw)-mf_fitted
          }
          y=residMF
          res=c(res,list(y));fittedFM=c(fittedFM,list(mf_fitted))
        }
      }
      mfTrainfd=fReg
    }
  }
  
  
  out=list(model=list(ml=ml,mf=temp),res=y,resList=res,'mfTrainfd'=mfTrainfd,fyl=fyList,'fittedFM'=fittedFM)
  return(out)
}


main3=function(response,lReg){
  ### this is for the case that the response are longitudinal data, while the covariates are
  ### scalers. 
  
  y=response
  ml=NULL
  if(!is.null(lReg)){
    ml=lm(y~lReg)
    y=as.matrix(resid(ml))
  }

  out=list(model=list(ml=ml,fl=NULL),res=y)
  return(out)
}




## creat fd object ##
mat2fd=function(mat,fdList=NULL){
  fl=list(time=seq(0,1,len=ncol(mat)),nbasis=min(as.integer(ncol(mat)/5),23),norder=6,bSpline=TRUE,Pen=c(0,0),lambda=1e-4)
  nbasis=c(fdList$nbasis,fl$nbasis)[1]
  norder=c(fdList$norder,fl$norder)[1]
  lambda=c(fdList$lambda,fl$norder)[1]
  bSpline=c(fdList$bSpline,fl$bSpline)[1]
  time=list(a=fdList$time,b=fl$time)
  time=time[[which(unlist(lapply(time,is.null))^2==0)[1]]]
  if(1-bSpline) fl$Pen=c(c(0,(2*pi/diff(range(time)))^2,0))
  Pen=list(a=fdList$Pen,b=fl$Pen)
  Pen=Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis=create.bspline.basis(range(time),nbasis,norder)
  if(1-bSpline) basis=create.fourier.basis(range(time),nbasis,diff(range(time)))
  Par=vec2Lfd(Pen,range(time))
  matfd=smooth.basisPar(time,t(mat),basis,Lfdobj=Par,lambda)$fd
  return(matfd)
}

betaPar=function(betaList=NULL){
  bl=list(rtime=c(0,1),nbasis=19,norder=4,bSpline=TRUE,Pen=c(0,0),lambda=1e4,bivar=FALSE,lambdas=1)
  nbasis=c(betaList$nbasis,bl$nbasis)[1]
  norder=c(betaList$norder,bl$norder)[1]
  lambda=c(betaList$lambda,bl$lambda)[1]
  bSpline=c(betaList$bSpline,bl$bSpline)[1]
  bivar=c(betaList$bivar,bl$bivar)[1]
  rtime=list(a=betaList$rtime,b=bl$rtime)
  rtime=rtime[[which(unlist(lapply(rtime,is.null))^2==0)[1]]]
  if(1-bSpline) fl$Pen=c(c(0,(2*pi/diff(rtime))^2,0))
  Pen=list(a=betaList$Pen,b=bl$Pen)
  Pen=Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis=create.bspline.basis(rtime,nbasis,norder)
  if(1-bSpline) basis=create.fourier.basis(rtime,nbasis,diff(rtime))
  Par=vec2Lfd(Pen,rtime)
  if(1-bivar) betaPar=fdPar(basis, Par, lambda)
  if(bivar){
    lambdas=c(betaList$lambdas,bl$lambdas)
    betaPar=bifdPar(bifd(matrix(0,nbasis,nbasis), basis, basis),
                    Par, Par, lambda, lambdas)
  } 
  return(betaPar)  
}


repgp.loglikelihood=function(hyper.p,response,Data,Cov,gamma=1,time=NULL,...){
  ### response is expected to be matrices with ncol replications and nrow observations
  ### Data is expected to be a list with matrices
  single=rep(0,ncol(response))
  for(i in 1:ncol(response)){
    old_input=as.matrix(do.call('cbind',lapply(Data,function(j) j=j[,i])))
    single[i]=gp.loglikelihood2(hyper.p,Data=old_input,response=response[,i],Cov=Cov,gamma=gamma,...)  
  }
  out=sum(single)
  return(out)
}




repgp.Dloglikelihood=function(hyper.p, response,Data,Cov,gamma=1,time=NULL,...){#,Xprior=prior_D1_likelihood,Xprior2=prior_likelihood){
  ### response is expected to be matrices with ncol replications and nrow observations
  ### Data is expected to be a list with matrices
  single=matrix(0,ncol=length(hyper.p),nrow=ncol(response))
  for(i in 1:nrow(single)){
    old_input=as.matrix(do.call('cbind',lapply(Data,function(j) j=j[,i])))
    single[i,]=gp.Dlikelihood2(hyper.p,Data=old_input,response=response[,i],Cov=Cov,gamma=gamma,...)  
  }
  out=apply(single,2,sum)
}


fisherinfo=function(pp.cg,X,Y,Cov,gamma){
  n=length(Cov)
  CovList=vector('list',n)
  for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(pp.cg,X,X,gamma=gamma))
    if(j!='cov.pow.ex')
      return(f(pp.cg,X,X))
  }  )
  if(length(CovL)==1)
    Q=CovL[[1]]
  if(length(CovL)>1)
    Q=Reduce('+',CovL)
  
  response=as.matrix(Y)
  X=as.matrix(X)
  Q=Q+diag(exp(pp.cg$vv),dim(Q)[1])
  invQ=mymatrix2(Q)$res
  QR=invQ%*%response
  AlphaQ=QR%*%t(QR)-invQ
  
  D2=function(d1,d2,inv.Q,Alpha.Q){
    Aii=t(d1)%*%inv.Q%*%d1
    al=Alpha.Q+inv.Q
    return(0.5*(sum(Alpha.Q*(d2-Aii))-sum(al*Aii)))
  }
  
  D2fx=lapply(seq_along(pp.cg),function(i){
    Dp=pp.cg[i]
    name.Dp=names(Dp)
    f=get(paste0('D2',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      D2para=f(pp.cg,X,gamma=gamma,inv.Q=invQ,Alpha.Q=AlphaQ)
    if(!name.Dp%in%c('pow.ex.w','pow.ex.v'))
      D2para=f(pp.cg,X,inv.Q=invQ,Alpha.Q=AlphaQ)
    return(D2para)
  })
  names(D2fx)=names(pp.cg)
  II=abs(-1/(unlist(D2fx)*dim(X)[1]))
  return(II)
}

gpfrtrain=function(response,lReg=NULL,fReg=NULL,fyList=NULL,fbetaList_l=NULL,fxList=NULL,
                   fbetaList=NULL,concurrent=TRUE,fbetaList_f=NULL,gpReg=NULL,hyper=NULL,Cov,gamma=1,
                   time=NULL,accuracy=c('high','normal','low'),trace.iter=5,fitting=FALSE,rPreIdx=FALSE){
  y=y_raw=response
  if(is.vector(y)) model='main1'
  if(is.matrix(y)){
    if(ncol(y)==1 | nrow(y)==1) model='main1'
    else model='main2'
  }
  if(is.fd(y) | is.fdata(y)) model ='main2'
  
  if(is.data.frame(y)) model='main3'
  ModelType=model
  
  if(model=='main1') 
    model=main1(response=response,lReg=lReg,fReg=fReg,fxList=fxList,fbetaList=fbetaList)
  if(model=='main2') 
    model=main2(response=response,lReg=lReg,fReg=fReg,fyList=fyList,fbetaList_l=fbetaList_l,
                fxList=fxList,concurrent=concurrent,fbetaList_f=fbetaList_f,time=time)  
  iuuL=NULL
  iuuL=mymatrix2(crossprod(cbind(lReg)))$res
  fittedFM=model$fittedFM
  
  ## convert fd/fdata class to matrix
  response=model$res;Data=gpReg
  if(class(response)!='matrix') stop('expecting matrix residual from functinoal regression')
  ftime=model$fyl$time
  if(!is.null(ftime))time=ftime
  if(is.null(ftime) & is.null(time)) stop('expecting input time')
  
  
  
  if(unique(unlist(lapply(Data,class)))=='fdata') Data=lapply(Data,function(i) i=t(i$data))

  if(class(response)=='matrix') response=t(response)
  if(unique(unlist(lapply(Data,class))=='matrix')) Data=lapply(Data,t)
  
  if(class(Data)=='fd')
    Data=(eval.fd(time,Data))
  if(unique(unlist(lapply(Data,class))=='fd')) Data=lapply(Data,function(i) (eval.fd(time,i)))
  
  
  ### this is an approximation of iuu for functional regression
  iuuF=NULL
  if(ModelType=='main2' & !is.null(fReg)){
    if(concurrent){
      iuuF=lapply(model$model$mf,function(i){
        BetaBasis=i$betaestlist$x$fd$basis
        xBasis=i$xfdlist$x$basis
        out=inprod(BetaBasis,xBasis)
        C=i$xfdlist$x$coefs
        out=out%*%C/(ncol(out)*nrow(C))
        out=tcrossprod(out)
        i=out
      })
    }
  }
 
  
  stepsize = nrow(response)  #training data size
  tdsize = as.integer(stepsize/2)  #choose half data for training
  
  
  sample_idx=1:tdsize*2-1
  response_raw=response; Data_raw=Data ## keep original data befreo reducing the dimension
  response=response[sample_idx,]
  if(ncol(Data[[1]])!=ncol(response)) Data=lapply(Data, t)
  Data=lapply(Data,function(i) i=i[sample_idx,])
  
  init0 = unlist(hyper)
  accuracy=accuracy[1]
  if(accuracy=='high') acc=1e-10
  if(accuracy=='normal') acc=1e-6
  if(accuracy=='low') acc=1e-2
  if(rPreIdx)   optm.idx=1:as.integer(nrow(response)/4)*3
  optm.idx=1:as.integer(nrow(response)/4)*3;
#   init1=nlminb(init0,repgp.loglikelihood,repgp.Dloglikelihood,response=response[optm.idx,],Data=lapply(Data,function(i) i=i[optm.idx,]),
#                Cov=Cov,gamma=gamma,control=list(iter.max=5,rel.tol=1e-1))[[1]]
  cat('    optimizing    ','\n')
#   hhhhp<<-init0;dadada<<-Data;rpp<<respo
  pp=nlminb(init0,repgp.loglikelihood,repgp.Dloglikelihood,response=response,Data=Data,Cov=Cov,gamma=gamma,
            control=list(trace=trace.iter,rel.tol=acc))#,Xprior=prior,Xprior2=NA)
  cat('    optimization done','\n','\n')
  pp.cg=pp[[1]]
  
  names(pp.cg)=names(init0)
  pp.df=data.frame(pp.cg=pp.cg,pp.N=substr(names(init0),1,8))
  names(pp.df)=c('pp.cg','pp.N')
  pp.cg=split(pp.df$pp.cg,pp.df$pp.N)
  
  allbat=vector('list',length=ncol(response))
  for(i in 1:ncol(response)){
    allbat[[i]]=cbind(response[,i],do.call('cbind',lapply(Data,function(j) j=j[,i])))
  }
  
  Qlist=lapply(allbat,function(l) fisherinfo(pp.cg=pp.cg,X=as.matrix(l[,-1]),Y=l[,1],Cov=Cov,gamma=gamma))
  II=abs(-1/apply(do.call('cbind',Qlist),1,sum))
  cat("fisher's information done",'\n','\n')
  
  mean=t(Reduce('+',fittedFM)) ## mean from functional regression model
  
  fitted=fitted.sd=NULL
  n=length(Cov);hyper.cg=pp.cg
  if(fitting){
    fitted=fitted.sd=matrix(0,ncol=ncol(response_raw),nrow=nrow(response_raw))
    for(i in seq_along(allbat)){
      dr=as.matrix(do.call('cbind',lapply(Data_raw,function(j) j=j[,i])))
      yy=as.matrix(response_raw[,i])
      
      CovList=vector('list',n)
      for(k in 1:n) CovList[k]=list(paste0('cov.',Cov[k]))
      
      CovL=lapply(CovList,function(j){
        f=get(j)
        if(j=='cov.pow.ex')
          return(f(hyper.cg,dr,dr,gamma=gamma))
        if(j!='cov.pow.ex')
          return(f(hyper.cg,dr,dr))
      }  )
      if(length(CovL)==1)
        Q=CovL[[1]]
      if(length(CovL)>1)
        Q=Reduce('+',CovL)
      
      Q=Q+diag(exp(hyper.cg$vv),dim(Q)[1])
      invQ=mymatrix2(Q)$res
      QR=invQ%*%yy
      AlphaQ=QR%*%t(QR)-invQ
      yfit=(Q-diag(exp(hyper.cg$vv),dim(Q)[1]))%*%invQ%*%(yy)+mean[,i]
      s2=exp(hyper.cg$vv)*rowSums((Q-diag(exp(hyper.cg$vv),dim(Q)[1]))*t(invQ)) ### ???
      fitted[,i]=yfit
      fitted.sd[,i]=sqrt(s2)
      if(i%%5==0) cat('fitting',i,' th curve','\n')
    }
    
  }
    
  result=list('hyper'=pp.cg,'I'=II, 'modellist'=model$model,'CovFun'=Cov,'gamma'=gamma,'init_resp'=y_raw,meanFM=mean,
              'resid_resp'=response_raw,'fitted.mean'=fitted,'fitted.sd'=fitted.sd,'ModelType'=ModelType,'lTrain'=lReg ,
              'fTrain'=fReg,'mfTrainfd'=model$mfTrainfd,'gpTrain'=Data_raw,'time'=time,'iuuL'=iuuL,'iuuF'=iuuF,
              'fittedFM'=fittedFM,'fyList'=fyList)
  class(result)='gpfr'
  return(result)
}

gpfrpred=function(object,TestData,NewTime=NULL,lReg=NULL,fReg=NULL,gpReg=NULL,GP_predict=TRUE){
  if(class(object)!='gpfr') stop('The object is expected to be a gpfda object','\n')
  
  model=object$modellist
  if(is.null(model$ml) & !is.null(lReg)){
    cat('    model with scaler variable is not found, ignoring lReg','\n')
    lReg=NULL
  }
  if(is.null(model$mf[[1]]) & !is.null(fReg)){
    cat('    model with functional variable is not found, ignoring lReg','\n')
    fReg=NULL
  }
  if(!is.null(model$ml) & is.null(lReg) & object$ModelType=='main1'){
    stop('    expecting input variable for model with scaler variable. ','\n')
    lReg=NULL
  }
  if(!is.null(model$ml) & is.null(lReg) & object$ModelType=='main2'){
    stop('    expecting input variable for model with scaler variable. ','\n')
    lReg=NULL
  }
  if(!is.null(model$mf[[1]]) & is.null(fReg)){
    stop('    expecting input variable for model with functional variable. ','\n')
    fReg=NULL
  }
  
  if(!is.null(model$ml) | !is.null(model$mf))
    rtime=c(model$ml$yhatfdobj$fd$basis$rangeval,model$mf[[1]]$yhatfdobj$basis$rangeval,
            model$mf[[1]]$yhatfdobj$fd$basis$rangeval)[1:2]
  if(is.null(model$ml) & is.null(model$mf)) rtime=c(0,1)
  
  if(class(TestData)=='matrix'){
    test=TestData
    test=t(test)
    if(is.null(NewTime)) time=seq(rtime[1],rtime[2],len=col(test))
    if(!is.null(NewTime)) time=NewTime
  }
  
  if(class(TestData)=='fd'){
    if(is.null(NewTime)) time=seq(rtime[1],rtime[2],len=TestData$fdnames$time)
    if(!is.null(NewTime)) time=NewTime
    test=eval.fd(time,TestData)
  } 
  testtime=time
  if(is.null(lReg)) lRegList=NULL
  if(!is.null(lReg)) lRegList=list(f=lReg)
  timeList=list(f=time)
  if(is.null(fReg)) fRegList=NULL
  if(!is.null(fReg))  fRegList=list(f=fReg)
  
  ml_var=0;mf_var=0
  
  if(!is.null(gpReg)& class(gpReg)!='list'){
    cat('Type I prediction is expecting gpReg to be a list with a response and an input. do type II prediction instead')
    gpReg=NULL
  }
  type=2
  if(!is.null(gpReg) & class(gpReg)=='list'){
    type=1
    nl=names(gpReg)
    if(sum(c('response','input','time')%in%names(gpReg))!=3)
      stop('check the name of the gpReg list. there must be one "response", one "input" and one "time"')
    if(!1%in%dim(as.matrix(gpReg$response)))
      stop('check the diemsion of the "response", it must be a column or row matrix, or a vector')
    gplReg=lReg
    gpfReg=gpReg$input
    gpresp=gpReg$response
    gptime=gpReg$time
    yhat_ml_list=vector('list',length=2)
    if(!is.null(lReg)) lRegList=list(f=lReg,gp=gplReg)
    if(!is.null(time)) timeList=list(f=time,gp=gptime)
    if(!is.null(fReg)) fRegList=list(f=fReg,gp=gpfReg)
    # else stop('expecting Test input')
  }
  
  ### predict ml for main2
  if(is.null(lRegList)) yhat_ml=0
  if(!is.null(lRegList)){
    for(ii in seq_along(lRegList)){
      lReg=lRegList[[ii]]
      time=timeList[[ii]]
      
      if (object$ModelType=='main2'){
        if(!is.null(lReg)){
          if(is.vector(lReg)) lReg=t(as.matrix(lReg))
          if(is.matrix(lReg)){
            if(length(model$ml$betaestlist)!=ncol(lReg))
              stop('dimension of lReg does not match the model')
            if(length(model$ml$betaestlist)==ncol(lReg)){
              x=model$ml$betaestlist
              for(i in 1:ncol(lReg))
                x[[i]]=as.matrix(lReg[,i])
            }
          }
          if(class(lReg)=='list'){
            if(unique(unlist(lapply(lReg,class)))=='fd')
              x=lapply(lReg,function(i) t(i$coefs))
            if(unique(unlist(lapply(lReg,class)))=='matrix')
              x=lReg
          }
          betalist=lapply(model$ml$betaestlist,function(i){
            i=predict(i,time)
          })
          if(ii==1) ml_var=do.call('cbind',x)%*%object$iuuL%*%t(do.call('cbind',x))
          for(i in 1:length(x)){
            x[[i]]=x[[i]]%*%t(betalist[[i]])
          }
          if(ii==1) yhat_ml=t(Reduce('+',x))
          if(ii==2) gpyhat_ml=t(Reduce('+',x))
        }
      }
    }
  }

  ### predict mf for main2
  if(is.null(fRegList)) yhat_mf=gpyhat_mf=matrix(0,ncol=1,nrow=1)
  if(!is.null(fRegList)){
    for(ii in seq_along(fRegList)){
      fReg=fRegList[[ii]]
      time=timeList[[ii]]
      if (object$ModelType=='main2'){
        if(!is.null(fReg)){
          if(is.fdata(fReg)) fReg=list(t(fReg$data))
          if(is.matrix(fReg) | is.fd(fReg)) fReg=list((fReg))
          if(unique(unlist(lapply(fReg,class)))=='fd'){
            fReg=lapply(fReg,function(i) eval.fd(time,i))
          }
          if(unique(unlist(lapply(fReg,class)))=='matrix'){
            fReg=lapply(fReg,t)
          }
          if(unique(unlist(lapply(fReg,ncol)))>1){
            if(length(fReg)!=1 & ii==1) stop('new samples of functional covaraites for functional regression are having wrong dimensions')
            if(length(fReg)!=1 & ii==2) stop('input functional covaraites for gaussian process are having wrong dimensions')
            if(length(fReg)==1){
              ftmp=vector('list',length=ncol(fReg[[1]]))
              for(i in seq_along(ftmp)) ftmp[[i]]=as.matrix(fReg[[1]][,i])
              fReg=ftmp
              rm(ftmp)
            }
          }
          
          fReg=lapply(fReg,function(i) i=cbind(matrix(1,nrow=nrow(fReg[[1]])),i) )
          ## find beta and multiply with the fx
          if(!'bifd'%in%unlist(lapply(model$mf[[1]],class))){
            fbeta=lapply(model$mf,function(i){
              intcept=predict(i$betaestlist[[1]],time)
              slope=predict(i$betaestlist[[2]],time)
              i=cbind(intcept,slope)
            })
            
            if(ii==1){
              mf_var=rep(0,length(fbeta))
              for(i in seq_along(fbeta)){
                iuuTime=seq(object$mfTrainfd[[i]]$basis$rangeval[1],object$mfTrainfd[[i]]$basis$rangeval[2],len=nrow(fReg[[i]]))
                Basis=object$mfTrainfd[[i]]$basis
                fx=as.matrix(fReg[[i]][,2])
                fx=t(eval.fd(seq(iuuTime[1],iuuTime[2],len=object$modellist$mf[[1]]$betaestlist$x$fd$basis$nbasis),
                             smooth.basis(iuuTime,fx,Basis)$fd))
                mf_var[i]=fx%*%object$iuuF[[i]]%*%t(fx)
              }
            }
            fReg[[i]]=apply(fReg[[i]],2,function(j){
              j=as.matrix(fbeta[[i]][,1])+as.matrix(fbeta[[i]][,2])*j
            })
            if(ii==1) yhat_mf=Reduce('+',fReg)
            if(ii==2) gpyhat_mf=Reduce('+',fReg)
          }
          
          if('bifd'%in%unlist(lapply(model$mf[[1]],class))){
            fbeta=lapply(model$mf,function(i){
              intcept=eval.fd(time,i$beta0estfd)
              slope=eval.bifd(time,time,i$beta1estbifd)
              i=cbind(intcept,slope)
            })
            for(i in seq_along(fbeta)){
              fReg[[i]]=apply(fReg[[i]],2,function(j){
                j=as.matrix(fbeta[[i]][,1])+as.matrix(fbeta[[i]][,-1])%*%as.matrix(j)
              })}
            if(i==1) yhat_mf=Reduce('+',fReg)
            if(i==2) gpyhat_mf=Reduce('+',fReg)
          }
          
        }
      }  
      
    }
  }
#   
  f.mean=yhat_ml+apply(yhat_mf,1,sum)
  f.var=apply((object$init_resp-t(object$resid_resp))^2,2,sum)/(nrow(object$init_resp)-1)
  object$fyList$time=object$time
  f.var=mat2fd(t(as.matrix(f.var)),object$fyList)
  f.var=abs(eval.fd(testtime,f.var))
  Fypredup = f.mean + 1.96*sqrt(f.var)
  Fypredlo = f.mean - 1.96*sqrt(f.var)
  if(1-GP_predict) return(list(ypred=cbind(f.mean,Fypredup,Fypredlo),unclass(object)))
  ## type I prediction
  
  
  if(type==1){
    
    cat('    Working out type I prediction','\n')
    gpresp=as.matrix(gpReg$response)
    if(nrow(gpresp)==1) gpresp=t(gpresp)
    ygpobs=as.matrix(gpresp)-gpyhat_ml-apply(gpyhat_mf,1,sum)    
    y_gppred=gppredict(train=FALSE,hyper=object$hyper, Data=as.matrix(gpReg$input), Y=as.matrix(ygpobs), Data.new=t(as.matrix(test)), Cov=object$CovFun,gamma=object$gamma)
    ygppred=y_gppred$pred.mean
    s2 = ((y_gppred$pred.sd)^2-exp(object$hyper$vv))%*%(1 + ml_var)#+sum(mf_var))
    ypred = yhat_ml+apply(yhat_mf,1,sum) + ygppred ## fd rgression plus gp regression
  }
  
  ## type II prediction
  if(is.null(gpReg)|type==2){
    cat('    Working out type II prediction','\n')
    
    fitted=matrix(0,ncol=ncol(object$resid_resp),nrow=ncol(test))
    fitted.var=fitted
    
    for(i in 1:ncol(object$resid_resp)){
      input=do.call('cbind',lapply(object$gpTrain,function(j) j=j[,i]))
      y_gppred=gppredict(train=FALSE,Data.new=t(as.matrix(test)),hyper=object$hyper,Data=as.matrix(input), Y=as.matrix(object$resid_resp[,i]),Cov=object$CovFun,gamma=object$gamma)
      ygppred = y_gppred$pred.mean
      s2 = ((y_gppred$pred.sd)^2-exp(object$hyper$vv))#%*%(1 + ml_var+sum(mf_var))
      fitted[,i]=yhat_ml+apply(yhat_mf,1,sum) + ygppred
      fitted.var[,i]=s2
    }
    ypred=as.matrix(apply(fitted,1,mean))
    s2=as.matrix(apply(fitted.var,1,mean))+apply(fitted^2,1,mean)-ypred^2
  }
  
  
  ypredup = ypred + 1.96*sqrt(s2)
  ypredlo = ypred - 1.96*sqrt(s2)
  
  CI=cbind(ypred, ypredup, ypredlo)
  
  result=c(list(ypred=CI, testtime=time,predtime=testtime,ypred.mean=ypred,ypred.sd=sqrt(s2)),unclass(object))
  class(result)='gpfr'
  return(result)
}



##### unfinished #####
gpfr=function(response,lReg=NULL,fReg=NULL,fyList=NULL,fbetaList_l=NULL,fxList=NULL,
              fbetaList=NULL,concurrent=TRUE,fbetaList_f=NULL,gpReg=NULL,hyper=NULL,Cov=c('pow.ex','linear'),gamma=1,
              time=NULL,NewHyper=NULL,accuracy=c('high','normal','low'),trace.iter=5,fitting=FALSE,rPreIdx=FALSE){
  if(is.list(gpReg)) col.no=length(gpReg)
  if(is.matrix(gpReg)) col.no=1
  if(!is.matrix(gpReg) & !is.list(gpReg)){
    cat('No gpReg found, doing functional regression only')
    col.no=1
  } 
  if(is.null(hyper)){
    hyper=list()
    if(any(Cov=='linear'))
      hyper$linear.a=rnorm(col.no)
    if(any(Cov=='pow.ex')){
      hyper$pow.ex.v=runif(1,-1,1)
      hyper$pow.ex.w=(-abs(rnorm(col.no)))
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.w=rnorm(col.no)
      hyper$rat.qu.s=runif(1,0.01,1)
      hyper$rat.qu.a=runif(1,0.01,1)
    }
    hyper$vv=sample(x=c(0.2,0.5),1)
    hyper.nam=names(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam=c(hyper.nam,NewHyper)
      nh.length=length(NewHyper)
      for(i in 1:nh.length){
        hyper=c(hyper,runif(1,-1,1))
      }
      names(hyper)=hyper.nam
    }
  }
  
  a1<-gpfrtrain(response=response,lReg=lReg,fReg=fReg,gpReg=gpReg,fyList=fyList,fbetaList_l=fbetaList_l,fxList=fxList,fbetaList_f=fbetaList_f,fbetaList=fbetaList,hyper=hyper,Cov=Cov,gamma=gamma,fitting=fitting,time=time,rPreIdx=rPreIdx,accuracy=accuracy,trace.iter=trace.iter,concurrent=concurrent)
  return(a1)
}

plot.gpfr=function (x, ..., type=c('raw','fitted','prediction')) 
{
  obj = x
  if(!type%in%c('raw','fitted','prediction')) 
    stop('type must be one of the raw, fitted or prediction')
  type=type[1]
  
  if(type=='raw'){
    matplot(obj$time,t(obj$init_resp),type='pl',pch=4,lwd=2,col='red',lty=3,main='Training',xlab='time')
  }
  
  if(type=='fitted'){
    matplot(obj$time,t(obj$init_resp),type='pl',pch=4,lwd=2,col='pink',lty=3,main='Training',xlab='time')
    for(i in 1:ncol(obj$fitted.mean)){
      polygon(c(obj$time, rev(obj$time)), c((obj$fitted.mean-obj$fitted.sd*1.96)[,i], rev((obj$fitted.mean+obj$fitted.sd*1.96)[,i])), 
              col = rgb(127,127,127,80, maxColorValue = 255), border = NA)
    }
    matpoints(obj$time,t(obj$init_resp),pch=4,cex=0.5,col='red',lty=3)
    matlines(obj$time,obj$fitted.mean,type='l',pch=4,lwd=0.5,col=4)
  }
  
  if(type=='prediction'){
    matplot(obj$time,t(obj$init_resp),type='l',pch=4,lwd=2,col='pink',lty=3,main='Training',xlab='time')
    matpoints(obj$time,t(obj$init_resp),pch=4,cex=0.5,col='red',lty=3)
    polygon(c(obj$predtime, rev(obj$predtime)), c(obj$ypred[,2], rev(obj$ypred[,3])), 
            col = rgb(127,127,127,100, maxColorValue = 255), border = NA)
    lines(obj$predtime,obj$ypred[,1],col=4,lwd=2)
    
  }
}


