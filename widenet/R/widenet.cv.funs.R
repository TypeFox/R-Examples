## slightly modified from glmnet's cvcompute function
## (same function is also internal in relaxnet package)

.relax.cvcompute <- function (mat, weights, foldid, nlams) 
{
    wisum = tapply(weights, foldid, sum)
    nfolds = max(foldid)
    outmat = matrix(NA, nfolds, ncol(mat))
    good = matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] = NA
    for (i in seq(nfolds)) {
        mati = mat[foldid == i, , drop = FALSE]
        wi = weights[foldid == i]
        outmat[i, ] = apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
        good[i, seq(nlams[i])] = 1
    }
    N = apply(good, 2, sum)
    list(cvraw = outmat, weights = wisum, N = N)
}

## adapted from cv.elnet from glmnet package

.widenet.cv.elnet <- function(outlist,
                              lambda,
                              x,
                              y,
                              weights,
                              offset,
                              foldid,
                              type.measure,
                              grouped,
                              order,
                              colsBinary,
                              screened.in.indices) {
  typenames=c(deviance="Mean-Squared Error",mse="Mean-Squared Error",mae="Mean Absolute Error")
  if(type.measure=="default")type.measure="mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE)){
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure="mse"
  }
     if(!is.null(offset))y=y-drop(offset)
     predmat=matrix(NA,length(y),length(lambda))
    nfolds=max(foldid)
    nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      fitobj$offset=FALSE

      ## expand basis if necessary

      x.to.use <- x[which, screened.in.indices[[i]], drop=FALSE]
      
      if (order == 2 || order == 3) {

        colsBinary.screened <- colsBinary[screened.in.indices[[i]]]

        numBinary.screened <- sum(colsBinary.screened == 2)
        numNotBinary.screened <- sum(colsBinary.screened == 3)

        x2.screened <- expand.to.order.2(x.to.use,
                                         colsBinary.screened,
                                         numBinary.screened,
                                         numNotBinary.screened)


        if(order == 2) {

          x.to.use <-
            cbind(x.to.use,
                  x2.screened)[, rownames(fitobj$beta)]

        } else  { ## order == 3

          x3.screened <- expand.to.order.3(x.to.use,
                                           x2.screened,
                                           colsBinary.screened,
                                           numBinary.screened,
                                           numNotBinary.screened)
          
          x.to.use <-
            cbind(x.to.use,
                  x2.screened,
                  x3.screened)[, rownames(fitobj$beta)]
        }
      }

      preds=predict(fitobj, x.to.use)

      nlami=length(outlist[[i]]$lambda)
       predmat[which,seq(nlami)]=preds
      nlams[i]=nlami
    }

  N=length(y) - apply(is.na(predmat),2,sum)
  cvraw=switch(type.measure,
    "mse"=(y-predmat)^2,
    "deviance"=(y-predmat)^2,
    "mae"=abs(y-predmat)
    )
   if( (length(y)/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }
 if(grouped){
   cvob=cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }

  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
}

## adapted from cv.elnet from glmnet and from .relax.cv.elnet from
## relaxnet

.relax.widenet.cv.elnet <- function(outlist,
                                    lambda,
                                    x,
                                    y,
                                    weights,
                                    offset,
                                    foldid,
                                    type.measure,
                                    grouped,
                                    order,
                                    colsBinary,
                                    screened.in.indices) {
  typenames=c(deviance="Mean-Squared Error",mse="Mean-Squared Error",mae="Mean Absolute Error")
  if(type.measure=="default")type.measure="mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE)){
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure="mse"
  }
     if(!is.null(offset))y=y-drop(offset)
     predmat=matrix(NA,length(y),length(lambda))
    nfolds=max(foldid)
    nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]

      if(inherits(fitobj, "relaxnet.intercept.only")) {

        predmat[which, ] <- fitobj
        ## just fill it in with intercept value stored in fitobj

        nlami <- length(lambda)
        
      } else {
        
        fitobj$offset=FALSE


        ## the beta mat, intercepts (a0) and lambda should already
        ## have been subset so that only the lower lambda values
        ## are processed (no need to repeat the values that were
        ## already done in the main model)

      ## expand basis if necessary

      x.to.use <- x[which, screened.in.indices[[i]], drop=FALSE]
      
      if (order == 2 || order == 3) {

        colsBinary.screened <- colsBinary[screened.in.indices[[i]]]

        numBinary.screened <- sum(colsBinary.screened == 2)
        numNotBinary.screened <- sum(colsBinary.screened == 3)

        x2.screened <- expand.to.order.2(x.to.use,
                                         colsBinary.screened,
                                         numBinary.screened,
                                         numNotBinary.screened)


        if(order == 2) {

          x.to.use <-
            cbind(x.to.use,
                  x2.screened)[, rownames(fitobj$beta)]

        } else  { ## order == 3

          x3.screened <- expand.to.order.3(x.to.use,
                                           x2.screened,
                                           colsBinary.screened,
                                           numBinary.screened,
                                           numNotBinary.screened)
          
          x.to.use <-
            cbind(x.to.use,
                  x2.screened,
                  x3.screened)[, rownames(fitobj$beta)]
        }
      }        

      if(order == 1) x.to.use <- x.to.use[, rownames(fitobj$beta), drop = FALSE]
      
     
      preds <- predict(fitobj,
                       x.to.use)
      nlami=length(fitobj$lambda)
      predmat[which,seq(nlami)]=preds
    }
      nlams[i]=nlami
    }

  N=length(y) - apply(is.na(predmat),2,sum)
  cvraw=switch(type.measure,
    "mse"=(y-predmat)^2,
    "deviance"=(y-predmat)^2,
    "mae"=abs(y-predmat)
    )
   if( (length(y)/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }
 if(grouped){
   cvob=.relax.cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }

  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
}

## adapted from cv.lognet from glmnet

.widenet.cv.lognet <- function(outlist,
                      lambda,
                      x,
                      y,
                      weights,
                      offset,
                      foldid,
                      type.measure,
                      grouped,
                      order,
                      colsBinary,
                      screened.in.indices) {
  typenames=c(mse="Mean-Squared Error",mae="Mean Absolute Error",deviance="Binomial Deviance",auc="AUC",class="Misclassification Error")
  if(type.measure=="default")type.measure="deviance"
  if(!match(type.measure,c("mse","mae","deviance","auc","class"),FALSE)){
    warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
    type.measure="deviance"
  }

###These are hard coded in the Fortran, so we do that here too
  prob_min=1e-5
  prob_max=1-prob_min
  ###Turn y into a matrix
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
    }
  N=nrow(y)
  nfolds=max(foldid)
  if( (N/nfolds <10)&&type.measure=="auc"){
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    type.measure="deviance"
  }
  if( (N/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }


  if(!is.null(offset)){
    is.offset=TRUE
    offset=drop(offset)
  }else is.offset=FALSE
  predmat=matrix(NA,nrow(y),length(lambda))
   nlams=double(nfolds)
  for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      if(is.offset)off_sub=offset[which]

      
      ## expand basis if necessary

      x.to.use <- x[which, screened.in.indices[[i]], drop=FALSE]
      
      if (order == 2 || order == 3) {

        colsBinary.screened <- colsBinary[screened.in.indices[[i]]]

        numBinary.screened <- sum(colsBinary.screened == 2)
        numNotBinary.screened <- sum(colsBinary.screened == 3)

        x2.screened <- expand.to.order.2(x.to.use,
                                         colsBinary.screened,
                                         numBinary.screened,
                                         numNotBinary.screened)


        if(order == 2) {

          x.to.use <-
            cbind(x.to.use,
                  x2.screened)[, rownames(fitobj$beta)]

        } else  { ## order == 3

          x3.screened <- expand.to.order.3(x.to.use,
                                           x2.screened,
                                           colsBinary.screened,
                                           numBinary.screened,
                                           numNotBinary.screened)
          
          x.to.use <-
            cbind(x.to.use,
                  x2.screened,
                  x3.screened)[, rownames(fitobj$beta)]
        }
      }

      
      preds=predict(fitobj,x.to.use, offset=off_sub,type="response")
      nlami=length(outlist[[i]]$lambda)
      predmat[which,seq(nlami)]=preds
      nlams[i]=nlami
    }
   ###If auc we behave differently
  if(type.measure=="auc"){
    cvraw=matrix(NA,nfolds,length(lambda))
    good=matrix(0,nfolds,length(lambda))
    for(i in seq(nfolds)){
      good[i,seq(nlams[i])]=1
      which=foldid==i
      for(j in seq(nlams[i])){
        cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
      }
    }
    N=apply(good,2,sum)
    weights=tapply(weights,foldid,sum)
  }
   else{
    ##extract weights and normalize to sum to 1
    ywt=apply(y,1,sum)
    y=y/ywt
    weights=weights*ywt

    N=nrow(y) - apply(is.na(predmat),2,sum)
    cvraw=switch(type.measure,
    "mse"=(y[,1]-(1-predmat))^2 +(y[,2]-predmat)^2,
    "mae"=abs(y[,1]-(1-predmat)) +abs(y[,2]-predmat),
    "deviance"= {
      predmat=pmin(pmax(predmat,prob_min),prob_max)
      lp=y[,1]*log(1-predmat)+y[,2]*log(predmat)
      ly=log(y)
      ly[y==0]=0
      ly=drop((y*ly)%*%c(1,1))
      2*(ly-lp)
    },
    "class"=y[,1]*(predmat>.5) +y[,2]*(predmat<=.5)
    )
 if(grouped){
   cvob=cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }
  }
   cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
}



## adapted from glmnet's cv.lognet and relaxnet's .relax.cv.lognet

.relax.widenet.cv.lognet <- function(outlist,
                                     lambda,
                                     x,
                                     y,
                                     weights,
                                     offset,
                                     foldid,
                                     type.measure,
                                     grouped,
                                     order,
                                     colsBinary,
                                     screened.in.indices){
  typenames=c(mse="Mean-Squared Error",mae="Mean Absolute Error",deviance="Binomial Deviance",auc="AUC",class="Misclassification Error")
  if(type.measure=="default")type.measure="deviance"
  if(!match(type.measure,c("mse","mae","deviance","auc","class"),FALSE)){
    warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
    type.measure="deviance"
  }

###These are hard coded in the Fortran, so we do that here too
  prob_min=1e-5
  prob_max=1-prob_min
  ###Turn y into a matrix
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }
  N=nrow(y)
  nfolds=max(foldid)
  if( (N/nfolds <10)&&type.measure=="auc"){
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    type.measure="deviance"
  }
  if( (N/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }


  if(!is.null(offset)){
    is.offset=TRUE
    offset=drop(offset)
  }else is.offset=FALSE
  predmat=matrix(NA,nrow(y),length(lambda))
  nlams=double(nfolds)
  for(i in seq(nfolds)){
    which=foldid==i
    fitobj=outlist[[i]]

    if(inherits(fitobj, "relaxnet.intercept.only")) {

      if(is.offset)
        stop("Internal Error:\n",
             "haven't dealt with offsets for binomial intercept models yet")
      
      predmat[which, ] <- 1 / (1 + exp(-fitobj))
      ## just fill it in with intercept value stored in fitobj
      ## check this again to make sure I'm doing this right
      
      nlami <- length(lambda)
        
    } else {

      ## the beta mat, intercepts (a0) and lambda should already
      ## have been subset so that only the lower lambda values
      ## are processed (no need to repeat the values that were
      ## already done in the main model)

      ## expand basis if necessary

      x.to.use <- x[which, screened.in.indices[[i]], drop=FALSE]
      
      if (order == 2 || order == 3) {

        colsBinary.screened <- colsBinary[screened.in.indices[[i]]]

        numBinary.screened <- sum(colsBinary.screened == 2)
        numNotBinary.screened <- sum(colsBinary.screened == 3)

        x2.screened <- expand.to.order.2(x.to.use,
                                         colsBinary.screened,
                                         numBinary.screened,
                                         numNotBinary.screened)


        if(order == 2) {

          x.to.use <-
            cbind(x.to.use,
                  x2.screened)[, rownames(fitobj$beta)]

        } else  { ## order == 3

          x3.screened <- expand.to.order.3(x.to.use,
                                           x2.screened,
                                           colsBinary.screened,
                                           numBinary.screened,
                                           numNotBinary.screened)
          
          x.to.use <-
            cbind(x.to.use,
                  x2.screened,
                  x3.screened)[, rownames(fitobj$beta)]
        }
      }        

      if(order == 1) x.to.use <- x.to.use[, rownames(fitobj$beta), drop = FALSE]
      

      

      if(is.offset) off_sub=offset[which]
      preds <- predict(fitobj,
                       x.to.use,
                       offset=off_sub, type="response")
      nlami=length(fitobj$lambda)
      predmat[which,seq(nlami)]=preds
    }
    
    nlams[i]=nlami
  }
  ##If auc we behave differently
  if(type.measure=="auc") {
    cvraw=matrix(NA,nfolds,length(lambda))
    good=matrix(0,nfolds,length(lambda))
    for(i in seq(nfolds)){
      good[i,seq(nlams[i])]=1
      which=foldid==i
      for(j in seq(nlams[i])){
        cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
      }
    }
    N=apply(good,2,sum)
    weights=tapply(weights,foldid,sum)

  } else {

    ##extract weights and normalize to sum to 1
    ywt=apply(y,1,sum)
    y=y/ywt
    weights=weights*ywt

    N=nrow(y) - apply(is.na(predmat),2,sum)
    cvraw=switch(type.measure,
      "mse"=(y[,1]-(1-predmat))^2 +(y[,2]-predmat)^2,
      "mae"=abs(y[,1]-(1-predmat)) +abs(y[,2]-predmat),
      "deviance"= {
        predmat=pmin(pmax(predmat,prob_min),prob_max)
        lp=y[,1]*log(1-predmat)+y[,2]*log(predmat)
        ly=log(y)
        ly[y==0]=0
        ly=drop((y*ly)%*%c(1,1))
        2*(ly-lp)},
      "class"=y[,1]*(predmat>.5) +y[,2]*(predmat<=.5)
      )
    if(grouped){
      cvob=.relax.cvcompute(cvraw,weights,foldid,nlams)
      cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
    }
  }
  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
}
