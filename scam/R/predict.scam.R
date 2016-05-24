
## clone of the predict.gam(mgcv) .....


predict.scam <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,
                       block.size=1000,newdata.guaranteed=FALSE,na.action=na.pass,...) 
{
# This function is used for predicting from a SCAM. object is a scam object, newdata a dataframe to
# be used in prediction......
#
# Type == "link"     - for linear predictor
#      == "response" - for fitted values
#      == "terms"    - for individual terms on scale of linear predictor 
#      == "iterms"   - exactly as "terms" except that se's include uncertainty about mean for unconstrained smooths
#      == "lpmatrix" - for matrix mapping parameters to l.p.
# Steps are:
#  1. Set newdata to object$model if no newdata supplied
#  2. split up newdata into manageable blocks if too large
#  3. Obtain parametric model matrix (safely!)
#  4. Work through smooths calling prediction.matrix constructors for each term
#  5. Work out required quantities
# 
# The splitting into blocks enables blocks of compiled code to be called efficiently
# using smooth class specific prediction matrix constructors, without having to 
# build up potentially enormous prediction matrices.
# if newdata.guaranteed == TRUE then the data.frame is assumed complete and
# ready to go, so that only factor levels are checked for sanity.
# 
# if `terms' is non null then it should be a list of terms to be returned 
# when type=="terms". 
# if `object' has an attribute `para.only' then only parametric terms of order
# 1 are returned for type=="terms": i.e. only what termplot can handle.
#
# if no new data is supplied then na.action does nothing, otherwise 
# if na.action == "na.pass" then NA predictors result in NA predictions (as lm
#                   or glm)
#              == "na.omit" or "na.exclude" then NA predictors result in
#                       dropping
# if GC is TRUE then gc() is called after each block is processed

  if (type!="link"&&type!="terms"&&type!="iterms"&&type!="response"&&type!="lpmatrix")  
  { warning("Unknown type, reset to terms.")
    type<-"terms"
  }
  if (!inherits(object,"scam")) stop("predict.scam can only be used to predict from scam objects")

#  if (ncol(attr(object$terms,"factors")) == 1)
#      { if (max(newdata) > max(object$data[,attr(b$terms,"term.labels")])) 
#           stop("predict.scam can only be used for data within #the range of observed values, please use extrapolate.scam #otherwise")  }

  ## to mimic behaviour of predict.lm, some resetting is required ...
  if (missing(newdata)) na.act <- object$na.action else
  { if (is.null(na.action)) na.act <- NULL 
    else {
      na.txt <- "na.pass"
      if (is.character(na.action))
      na.txt <- substitute(na.action) else
      if (is.function(na.action)) na.txt <- deparse(substitute(na.action))
      if (na.txt=="na.pass") na.act <- "na.exclude" else
      if (na.txt=="na.exclude") na.act <- "na.omit" else
      na.act <- na.action
    }
  } ## ... done
  # get data from which to predict.....  
  nd.is.mf <- FALSE # need to flag if supplied newdata is already a model frame
  if (newdata.guaranteed==FALSE)
  { if (missing(newdata)) # then "fake" an object suitable for prediction 
    { newdata<-object$model
      new.data.ok <- FALSE
      nd.is.mf <- TRUE
    }
    else  # do an R ``standard'' evaluation to pick up data
    { new.data.ok <- TRUE
      if (is.data.frame(newdata)&&!is.null(attr(newdata,"terms"))) # it's a model frame
      { if (sum(!(names(object$model)%in%names(newdata)))) stop(
        "newdata is a model.frame: it should contain all required variables\n")
         nd.is.mf <- TRUE
      } else
      { ## Following is non-standard to allow convenient splitting into blocks
        ## below, and to allow checking that all variables are in newdata ...

        ## get names of required variables, less response, but including offset variable
        Terms <- delete.response(terms(object))
        allNames <- all.vars(Terms)
        ff <- reformulate(allNames)
        if (sum(!(allNames%in%names(newdata)))) { 
        warning("not all required variables have been supplied in  newdata!\n")}
        ## note that `xlev' argument not used here, otherwise `as.factor' in 
        ## formula can cause a problem ... levels reset later.
        newdata <- eval(model.frame(ff,data=newdata,na.action=na.act),parent.frame()) 
        na.act <- attr(newdata,"na.action")
      }
    }
  } else { ## newdata.guaranteed == TRUE
    na.act <- NULL
    new.data.ok=TRUE ## it's guaranteed!
  }
  

  if (new.data.ok)
  { ## check factor levels are right ...
    names(newdata)->nn # new data names
    colnames(object$model)->mn # original names
    for (i in 1:length(newdata)) 
    if (nn[i]%in%mn && is.factor(object$model[,nn[i]])) # then so should newdata[[i]] be 
    { newdata[[i]]<-factor(newdata[[i]],levels=levels(object$model[,nn[i]])) # set prediction levels to fit levels
    }

    # split prediction into blocks, to avoid running out of memory
    if (length(newdata)==1) newdata[[2]]<-newdata[[1]] # avoids data frame losing its labels and dimensions below!
    if (is.null(dim(newdata[[1]]))) np<-length(newdata[[1]]) 
    else np<-dim(newdata[[1]])[1] 
    nb<-length(object$coefficients)
    if (block.size<1) block.size <- np
    n.blocks<-np%/%block.size
    b.size<-rep(block.size,n.blocks)
    last.block<-np-sum(b.size)
    if (last.block>0) 
    { n.blocks<-n.blocks+1  
      b.size[n.blocks]<-last.block
    }
  } else # no new data, just use object$model
  { np <- nrow(object$model)
    nb <- length(object$coefficients)
    n.blocks <- 1
    b.size <- array(np,1)
  }
  # setup prediction arrays
  n.smooth<-length(object$smooth)
  if (type=="lpmatrix")
  { H<-matrix(0,np,nb)
  } else
  if (type=="terms"||type=="iterms")
  { term.labels<-attr(object$pterms,"term.labels")
    if (is.null(attr(object,"para.only"))) para.only <-FALSE else
    para.only <- TRUE  # if true then only return information on parametric part
    n.pterms <- length(term.labels)
    fit<-array(0,c(np,n.pterms+as.numeric(!para.only)*n.smooth))
    if (se.fit) se<-fit
    ColNames<-term.labels
  } else
  { fit<-array(0,np)
    if (se.fit) se<-fit
  }
  stop<-0

  Terms <- delete.response(object$pterms)
  s.offset <- NULL # to accumulate any smooth term specific offset
  any.soff <- FALSE # indicator of term specific offset existence
  for (b in 1:n.blocks)  # work through prediction blocks
  { start<-stop+1
    stop<-start+b.size[b]-1
    if (n.blocks==1) data <- newdata else data<-newdata[start:stop,]
    X <- matrix(0,b.size[b],nb)
    Xoff <- matrix(0,b.size[b],n.smooth) ## term specific offsets 
    ## implements safe prediction for parametric part as described in
    ## http://developer.r-project.org/model-fitting-functions.txt
    if (new.data.ok)
    { if (nd.is.mf) mf <- model.frame(data,xlev=object$xlevels) else
      { mf <- model.frame(Terms,data,xlev=object$xlevels)
        if (!is.null(cl <- attr(object$pterms,"dataClasses"))) .checkMFClasses(cl,mf)
      } 
      Xp <- model.matrix(Terms,mf,contrasts=object$contrasts) 
    } else 
    { Xp <- model.matrix(Terms,object$model)
      mf <- newdata # needed in case of offset, below
    }
    
    if (object$nsdf) X[,1:object$nsdf]<-Xp
    if (n.smooth) for (k in 1:n.smooth) 
    { Xfrag <- PredictMat(object$smooth[[k]],data)	

### added code specific for scam....
	if (inherits(object$smooth[[k]], c("mpi.smooth","mpd.smooth", "cv.smooth", "cx.smooth",   "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth","tedmi.smooth","tedmd.smooth")))
           X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag[,2:ncol(Xfrag)]
      else if (inherits(object$smooth[[k]], c("tesmi1.smooth","tesmi2.smooth","tesmd1.smooth",
                  "tesmd2.smooth")))  # for single monotonicity...
             X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag%*%object$smooth[[k]]$Zc
      else
          X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
      Xfrag.off <- attr(Xfrag,"offset") ## any term specific offsets?
      if (!is.null(Xfrag.off)) { Xoff[,k] <- Xfrag.off; any.soff <- TRUE }
      if (type=="terms"||type=="iterms") ColNames[n.pterms+k]<-object$smooth[[k]]$label
    }
    # have prediction matrix for this block, now do something with it
    if (type=="lpmatrix") { 
      H[start:stop,]<-X
      if (any.soff) s.offset <- rbind(s.offset,Xoff)
    } else 
    if (type=="terms" ||type=="iterms")
    {
      ind <- 1:length(object$assign)
      if (n.pterms)  # work through parametric part
      for (i in 1:n.pterms)
      { ii <- ind[object$assign==i]

#### CORRECTIONS FOR SCAM....
        fit[start:stop,i] <- as.matrix(X[,ii])%*%object$coefficients.t[ii]
        if (se.fit) se[start:stop,i]<-
        sqrt(rowSums((as.matrix(X[,ii])%*%object$Vp.t[ii,ii])*as.matrix(X[,ii])))
      }

      if (n.smooth&&!para.only) 
      { for (k in 1:n.smooth) # work through the smooth terms 
        { first<-object$smooth[[k]]$first.para;last<-object$smooth[[k]]$last.para

### CORRECTED for SCAM ...

          fit[start:stop,n.pterms+k]<-X[,first:last]%*%object$coefficients.t[first:last] + Xoff[,k]
          if (se.fit) { # diag(Z%*%V%*%t(Z))^0.5; Z=X[,first:last]; V is sub-matrix of Vp
               if (inherits(object$smooth[[k]], c("mpi.smooth","mpd.smooth","cv.smooth", "cx.smooth", "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth"))){
                          if (nrow(X)==1) # prediction vector if prediction is made for only one value of covariates
                              X1 <- c(1,t(X[,first:last]))
                          else 
                              X1 <- cbind(rep(1,nrow(X)),X[,first:last]) # prediction matrix
                                # X0 - model matrix of the original data....
                          object.X <- model.matrix(object)
                          X0 <- cbind(rep(1,nrow(object.X)),object.X[,first:last]) 
                          q <- ncol(X0)
                          onet <- matrix(rep(1,nrow(X0)),1,nrow(X0))
                          A <- onet%*%X0
                          qrX <- qr(X0)
                          R <- qr.R(qrX) 
                          qrA <- qr(t(A))
                          R <- R[-1,]
                          RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                          RZa.inv <- solve(RZa)
                          RZaR <- RZa.inv%*%R
                          if (nrow(X)==1)
                              XZa <- t(qr.qty(qrA,X1))[,2:q]
                          else
                              XZa <- t(qr.qty(qrA,t(X1)))[,2:q]
                          Ga <- XZa%*%RZaR
                          Vp <- object$Vp.t[c(1,first:last),c(1,first:last)] 
                          Vp.c <- Vp
                          Vp.c[,1] <- rep(0,nrow(Vp))
                          Vp.c[1,] <- rep(0,ncol(Vp))
                          se[start:stop,n.pterms+k] <- sqrt(rowSums((Ga%*%Vp.c)*Ga))
               }
               else if (inherits(object$smooth[[k]], c("tedmi.smooth","tedmd.smooth", 
                         "tesmi1.smooth","tesmi2.smooth","tesmd1.smooth","tesmd2.smooth"))) {
                                # X0 - model matrix of the original data....
                        object.X <- model.matrix(object)
                        X0 <- cbind(rep(1,nrow(object.X)),object.X[,first:last]) 
                        onet <- matrix(rep(1,nrow(X0)),1,nrow(X0))
                        if (nrow(X)==1) # prediction vector if prediction is made for only one value of covariates
                              X1 <- c(1,t(X[,first:last]))
                          else 
                              X1 <- cbind(rep(1,nrow(X)),X[,first:last]) # prediction matrix
                        A <- onet%*%X0
                        qrX <- qr(X0)
                        R <- qr.R(qrX)
                        qrA <- qr(t(A))
                        q <- ncol(X0)
                       if (inherits(object$smooth[[k]], c("tedmi.smooth","tedmd.smooth")))
                          { # get RZaR for double monotonicity...
                            R <- R[-1,]
                            RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                            RZa.inv <- solve(RZa)
                            RZaR <- RZa.inv%*%R
                          }
                       else # get RZaR for single monotonicity...
                          { RZa <- t(qr.qty(qrA,t(R)))[,2:q]
                            RZatRZa.inv <- solve(crossprod(RZa)) ## solve(t(RZa)%*%RZa) 
                            Q <- qr.Q(qrX)
                            B1 <-  tcrossprod(RZatRZa.inv,RZa) ## RZatRZa.inv%*%t(RZa)
                            RZaR <- B1%*%R
                           }
                       if (nrow(X)==1)
                              XZa <- t(qr.qty(qrA,X1))[,2:ncol(X1)]
                       else
                              XZa <- t(qr.qty(qrA,t(X1)))[,2:ncol(X1)]
                       Ga <- XZa%*%RZaR%*%object$smooth[[k]]$Zc
                       Vp <- object$Vp.t[first:last,first:last]
                       se[start:stop,n.pterms+k] <- sqrt(rowSums((Ga%*%Vp)*Ga))
               }
               else { ## for unconstrained smooth terms..... 
                    if (type=="iterms"&& attr(object$smooth[[k]],"nCons")>0) { ## termwise se to "carry the intercept
                      X1 <- matrix(object$cmX,nrow(X),ncol(X),byrow=TRUE)
                      meanL1 <- object$smooth[[k]]$meanL1
                      if (!is.null(meanL1)) X1 <- X1 / meanL1              
                      X1[,first:last] <- X[,first:last]
                      se[start:stop,n.pterms+k] <- sqrt(pmax(0,rowSums((X1%*%object$Vp)*X1)))
                    } else  se[start:stop,n.pterms+k] <- ## terms strictly centred
                    sqrt(rowSums((X[,first:last]%*%object$Vp.t[first:last,first:last])*X[,first:last]))
               }         
          } ## end if (se.fit)
        }
        colnames(fit) <- ColNames
        if (se.fit) colnames(se) <- ColNames
      } else { # para.only
        colnames(fit) <- term.labels
        if (se.fit) colnames(se) <- term.labels
        # retain only terms of order 1 - this is to make termplot work
        order <- attr(object$pterms,"order")
        term.labels <- term.labels[order==1]
        fit <- as.matrix(as.matrix(fit)[,order==1])
        colnames(fit) <- term.labels
        if (se.fit) { se <- as.matrix(as.matrix(se)[,order==1])
        colnames(se) <- term.labels } 
      }
      if (!is.null(terms)) # return only terms requested via `terms'
      { if (sum(!(terms %in%colnames(fit)))) 
        warning("non-existent terms requested - ignoring")
        else { names(term.labels) <- term.labels
          term.labels <- term.labels[terms]  # names lost if only one col
          fit <- as.matrix(as.matrix(fit)[,terms])
          colnames(fit) <- term.labels
          if (se.fit) {se <- as.matrix(as.matrix(se)[,terms])
          colnames(se) <- term.labels}
        }
      }
    } else # "link" or "response"
    { k<-attr(attr(object$model,"terms"),"offset")

## CORRECTED for SCAM...
      fit[start:stop]<-X%*%object$coefficients.t + rowSums(Xoff)
      if (!is.null(k)) fit[start:stop]<-fit[start:stop]+model.offset(mf) + rowSums(Xoff)
      if (se.fit) se[start:stop]<-sqrt(rowSums((X%*%object$Vp.t)*X))
      if (type=="response") # transform    
      { fam<-object$family;linkinv<-fam$linkinv;dmu.deta<-fam$mu.eta  
        if (se.fit) se[start:stop]<-se[start:stop]*abs(dmu.deta(fit[start:stop])) 
        fit[start:stop]<-linkinv(fit[start:stop])
      }
    }
    rm(X)
   
  } ## end of prediction block loop
  rn <- rownames(newdata)
  if (type=="lpmatrix") { 
    colnames(H) <- names(object$coefficients);rownames(H)<-rn
    if (!is.null(s.offset)) { 
      s.offset <- napredict(na.act,s.offset)
      attr(H,"offset") <- s.offset
    }
    H <- napredict(na.act,H)
  } else { 
    if (se.fit) { 
      if (is.null(nrow(fit))) {
        names(fit) <- rn
        names(se) <- rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se) 
      } else { 
        rownames(fit)<-rn
        rownames(se)<-rn
        fit <- napredict(na.act,fit)
        se <- napredict(na.act,se)
      }
      H<-list(fit=fit,se.fit=se) 
    } else { 
      H<-fit
      if (is.null(nrow(H))) names(H) <- rn else
      rownames(H)<-rn
      H <- napredict(na.act,H)
    }
  }
  if (type=="terms"||type=="iterms") attr(H,"constant") <- object$coefficients[1]
  H # ... and return
}



###############################################################################
### ISSUES.....
#
#* predict.scam "terms" and "iterms" don't deal properly with 
#  shape constrained smooths: se.fit retured are the same for both types

