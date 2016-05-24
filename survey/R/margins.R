
marginpred<-function(model, adjustfor, predictat, ...) UseMethod("marginpred", model)

##
## Basic strategy: calibrate on ~model*adjustfor, to set interactions to zero
##

marginpred.svycoxph<-function(model, adjustfor, predictat, se=FALSE, ...){

  if(NROW(predictat)==0) return(NULL)
  
  design<-model$survey.design
  ##if (inherits(design,"twophase")) stop("Two-phase designs not yet supported")
  if (!is.null(model$na.action)) design<-design[-model$na.action,]

  modelformula<-formula(model)
  calformula<-eval(bquote( ~(.(modelformula[[3]]))*(.(adjustfor))))

  adjmf<-model.frame(terms(adjustfor), model.frame(design))
  adjmm<-model.matrix(terms(adjustfor), adjmf)
  modelmm<-model.matrix(model)[,-1,drop=FALSE]
  modelmm <- sweep(modelmm,2,model$means)

  if (qr(modelmm)$rank<ncol(modelmm))
    stop("model is singular")
  
  qrmain<-qr(cbind(modelmm,adjmm))
  if(qrmain$rank<(ncol(modelmm)+ncol(adjmm))){
    if (qrmain$rank==ncol(modelmm))
      stop("adjustment variables are all in model")
    adjmm<-adjmm[, (qrmain$pivot[-(1:ncol(modelmm))]-ncol(modelmm))[1:(qrmain$rank-ncol(modelmm))],drop=FALSE]
  }

  mm<-matrix(ncol=ncol(modelmm)*ncol(adjmm),nrow=nrow(modelmm))
  for(i in 1:ncol(modelmm)){
    mm[,(i-1)*ncol(adjmm)+(1:ncol(adjmm))]<-adjmm*modelmm[,i]
  }

  pop<-as.vector(outer(colSums(adjmm*weights(design)),
                       colSums(modelmm*weights(design))/sum(weights(design))
                       )
                 )

  g<-grake(mm,weights(design), calfun=cal.raking, bounds=c(0,Inf),
                population=pop, epsilon=1e-4,maxit=100,verbose=FALSE)
  
  if ( !is.null(attr(g,"failed"))) stop("Calibration failed")
  design$prob<-design$prob/g
  
  whalf<-sqrt(weights(design))
  tqr<-qr(mm*whalf)
  caldata <- list(qr=tqr, w=g*whalf, stage=0, index=NULL)
  
  class(caldata) <- c("greg_calibration","gen_raking")
  
  design$postStrata <- c(design$postStrata, list(caldata))
  design$call <- sys.call(-1)
  
  model<-eval(bquote(svycoxph(.(formula(model)), design=design)))
  predict(model, newdata=predictat, se=se,type="curve",...)

}


marginpred.svykmlist<-function(model, adjustfor, predictat, se=FALSE, ...){
  
  design<-eval.parent(attr(model, "call")$design)
  formula<-formula(model)
  if(!is.null(drop<-attr(model,"na.action")))
    design<-design[-drop,]
  
  if(NROW(predictat)==0) return(NULL)
  
  ##if (inherits(design,"twophase")) stop("Two-phase designs not yet supported")

  modelformula<-formula(model)
  calformula<-eval(bquote( ~(.(modelformula[[3]]))*(.(adjustfor))))

  adjmf <- model.frame(terms(adjustfor), model.frame(design),na.action=na.fail)
  adjmm <- model.matrix(terms(adjustfor), adjmf)
  adjmm <- sweep(adjmm, 2, colSums(adjmm*weights(design))/sum(weights(design)))
  modelmf <- model.frame(terms(formula), model.frame(design),na.action=na.fail)
  modelmm <- model.matrix(terms(formula), modelmf)[,-1,drop=FALSE]
  modelmm <- sweep(modelmm, 2, colSums(modelmm*weights(design))/sum(weights(design)))
  
  if (qr(modelmm)$rank<ncol(modelmm))
    stop("model is singular")
  
  qrmain<-qr(cbind(modelmm, adjmm))
  if(qrmain$rank<(ncol(modelmm)+ncol(adjmm))){
    if (qrmain$rank==ncol(modelmm))
      stop("adjustment variables are all in model")
    adjmm<-adjmm[, (qrmain$pivot[-(1:ncol(modelmm))]-ncol(modelmm))[1:(qrmain$rank-ncol(modelmm))],drop=FALSE]
  }

  mm<-matrix(ncol=ncol(modelmm)*ncol(adjmm),nrow=nrow(modelmm))
  for(i in 1:ncol(modelmm)){
    mm[,(i-1)*ncol(adjmm)+(1:ncol(adjmm))]<-adjmm*modelmm[,i]
  }

  pop<-as.vector(outer(colSums(adjmm*weights(design)),
                       colSums(modelmm*weights(design))/sum(weights(design))
                       )
                 )

  g<-grake(mm,weights(design), calfun=cal.raking, bounds=c(0, Inf),
                population=pop, epsilon=1e-4, maxit=100, verbose=FALSE)
  
  if ( !is.null(attr(g,"failed"))) stop("Calibration failed")
  design$prob<-design$prob/g
  
  whalf<-sqrt(weights(design))
  tqr<-qr(mm*whalf)
  caldata <- list(qr=tqr, w=g*whalf, stage=0, index=NULL)
  
  class(caldata) <- c("greg_calibration","gen_raking")
  
  design$postStrata <- c(design$postStrata, list(caldata))
  design$call <- sys.call(-1)
  
  eval(bquote(svykm(.(formula), design=design,se=.(se))))
}


marginpred.svyglm<-function(model, adjustfor, predictat,  ...){
  
  design<-model$survey.design
  formula<-formula(model)
  if(!is.null(drop<-attr(model,"na.action")))
    design<-design[-drop,]
  
  if(NROW(predictat)==0) return(NULL)
  
  ##if (inherits(design,"twophase")) stop("Two-phase designs not yet supported")

  modelformula<-formula(model)
  calformula<-eval(bquote( ~(.(modelformula[[3]]))*(.(adjustfor))))

  adjmf <- model.frame(terms(adjustfor), model.frame(design),na.action=na.fail)
  adjmm <- model.matrix(terms(adjustfor), adjmf)
  adjmm <- sweep(adjmm, 2, colSums(adjmm*weights(design))/sum(weights(design)))
  modelmm <- model.matrix(model)[,-1,drop=FALSE]
  modelmm <- sweep(modelmm, 2, colSums(modelmm*weights(design))/sum(weights(design)))
  
  if (qr(modelmm)$rank<ncol(modelmm))
    stop("model is singular")
  
  qrmain<-qr(cbind(modelmm, adjmm))
  if(qrmain$rank<(ncol(modelmm)+ncol(adjmm))){
    if (qrmain$rank==ncol(modelmm))
      stop("adjustment variables are all in model")
    adjmm<-adjmm[, (qrmain$pivot[-(1:ncol(modelmm))]-ncol(modelmm))[1:(qrmain$rank-ncol(modelmm))],drop=FALSE]
  }

  mm<-matrix(ncol=ncol(modelmm)*ncol(adjmm),nrow=nrow(modelmm))
  for(i in 1:ncol(modelmm)){
    mm[,(i-1)*ncol(adjmm)+(1:ncol(adjmm))]<-adjmm*modelmm[,i]
  }

  pop<-as.vector(outer(colSums(adjmm*weights(design)),
                       colSums(modelmm*weights(design))/sum(weights(design))
                       )
                 )

  g<-grake(mm,weights(design), calfun=cal.raking, bounds=c(0,Inf),
                population=pop, epsilon=1e-4, maxit=100, verbose=FALSE)
  
  if ( !is.null(attr(g,"failed"))) stop("Calibration failed")
  .design<-design
  .design$prob<-.design$prob/g
  
  whalf<-sqrt(weights(.design))
  tqr<-qr(mm*whalf)
  caldata <- list(qr=tqr, w=g*whalf, stage=0, index=NULL)
  
  class(caldata) <- c("greg_calibration","gen_raking")
  
  design$postStrata <- c(design$postStrata, list(caldata))
  design$call <- sys.call(-1)

  call<-model$call
  call$design<-quote(.design)
  newmodel<-eval(call)
  predict(newmodel,newdata=predictat, ...)
}

