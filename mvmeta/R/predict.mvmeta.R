###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
predict.mvmeta <-
function(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE,
  interval=c("confidence","prediction"), ci.level=0.95,
  format=c("matrix","list"), aggregate=c("stat","y"), na.action=na.pass, ...) {
#
################################################################################
# CREATE DESIGN MATRIX X AND OFFSET
#
  if(missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    offset <- object$offset
    # IN CASE, RE-INSERT MISSING VALUES FROM FITTING (ORIGINAL OBJECT)
    if(!is.null(object$na.action)) {
      X <- napredict(object$na.action,X)
      offset <- napredict(object$na.action,offset)
    }
  }  else {
    ttnr <- delete.response(terms(object))
    # HERE MISSING VALUES ARE HANDLES BY THE ARGUMENT OF predict
    mf <- model.frame(ttnr,newdata,na.action=na.action,xlev=object$xlevels)
    if(!is.null(class <- attr(ttnr,"dataClasses"))) .checkMFClasses(class,mf)
    X <- model.matrix(ttnr,mf,contrasts.arg=object$contrasts)
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset)!=NROW(X)) stop("wrong offset")
  }
#
################################################################################
# RE-CREATE OBJECTS
#
  interval <- match.arg(interval,c("confidence","prediction"))
  format <- match.arg(format,c("matrix","list"))
  aggregate <- match.arg(aggregate,c("stat","y"))
#
  # MISSING PATTERN (ONLY RELEVANT FOR PREDICTION ON ORIGINAL DATA)
  if(missing(newdata) || is.null(newdata)) {
    res <- as.matrix(object$residuals)
    if(!is.null(object$na.action)) res <- napredict(object$na.action,res)
  } else res <- matrix(0,nrow(X),object$dim$k)
  nay <- is.na(res)
#
  # TRANSFORM X IN LIST, ACCOUNTING FOR MISSING
  Xlist <- lapply(seq(nrow(X)),function(i) {
    dd <- diag(1,object$dim$k)
    diag(dd)[nay[i,]] <- NA
    return(dd%x%X[i,,drop=FALSE])})
#
################################################################################
# COMPUTE PREDICTION, SE, CI AND VCOV
#
  # COMPUTE PREDICTION, ACCOUNTING FOR OFFSET
  predlist <- lapply(seq(nrow(X)),function(i) {
    fit <- drop(Xlist[[i]]%*%as.numeric(object$coefficients))
    if(!is.null(offset)) fit <- fit + offset[i]
    names(fit) <- object$lab$k
    return(fit)})
#
  # COMPUTE VCOV, SE AND CI
  zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
  vcovlist <- lapply(Xlist,function(X) {
    vcov <- X %*% tcrossprod(object$vcov,X)
    if(interval=="prediction"&& object$method!="fixed") vcov <- vcov+object$Psi 
    dimnames(vcov) <- list(object$lab$k,object$lab$k)
    return(vcov)})
  selist <- lapply(vcovlist, function(x) sqrt(diag(x)))
  cilblist <- mapply(function(fit,se) fit-zvalci*se,
      predlist,selist,SIMPLIFY=FALSE)
  ciublist <- mapply(function(fit,se) fit+zvalci*se,
      predlist,selist,SIMPLIFY=FALSE)
#
################################################################################
# AGGREGATE AND RETURN
#
  # IF VCOV WHEN MORE THAN 1 OUTCOME, SWITCH TO LIST
  if(vcov && object$dim$k>1) format <- "list"
  # IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
  if(format=="matrix" && object$dim$k==1) {
    fit <- drop(do.call("rbind",predlist))
    names(fit) <- rownames(X)
    if(se) fit <- cbind(fit,se=drop(do.call("rbind",selist)))
    if(ci) fit <- cbind(fit,ci.lb=drop(do.call("rbind",cilblist)),
      ci.ub=drop(do.call("rbind",ciublist)))
    if(vcov) fit <- cbind(fit,vcov=drop(do.call("rbind",vcovlist)))
  # IF NO OTHER STAT, FORCE THE AGGREGATION ON PREDICTION
  } else if(format=="matrix" && !se && !ci) {
    fit <- do.call("rbind",predlist)
    rownames(fit) <- rownames(X)
    # WHEN AGGREGATE ON STAT
  } else if (format=="matrix" && aggregate=="stat") {
    fit <- list(fit=do.call("rbind",predlist))
    rownames(fit$fit) <- rownames(X)
    if(se) {
      fit$se <- do.call("rbind",selist)
      rownames(fit$se) <- rownames(X)
    }
    if(ci) {
      fit$ci.lb <- do.call("rbind",cilblist)
      fit$ci.ub <- do.call("rbind",ciublist)
      rownames(fit$ci.lb) <- rownames(fit$ci.ub) <- rownames(X)
    }
  # WHEN AGGREGATE ON OUTCOME
  } else if (format=="matrix" && aggregate=="y") {
    fit <- list()
    for(j in seq(object$dim$k)) {
      templist <- lapply(seq(predlist),function(i) {
        fit <- predlist[[i]][j]
        if(se) fit <- cbind(fit,selist[[i]][j])
        if(ci) fit <- cbind(fit,cilblist[[i]][j],
          ciublist[[i]][j])
        return(as.matrix(fit))})
      temp <- do.call("rbind",templist)
      dimnames(temp) <- list(rownames(X),
        c("fit","se","ci.lb","ci.ub")[c(TRUE,se,ci,ci)])
      fit[[j]] <- temp
    }
    names(fit) <- object$lab$k
  # ALL THE OTHER COMBINATIONS, PRODUCE A LIST
  } else {
    fit <- mapply(function(predarg,searg,cilbarg,ciubarg,vcovarg) {
      temp <- list(fit=predarg)
      if(se) temp$se <- searg
      if(ci) temp$ci.lb <- cilbarg
      if(ci) temp$ci.ub <- ciubarg
      if(vcov) temp$vcov <- vcovarg
      # RETURN ONE OBJECT IF ONLY FIT 
      if(se||ci||vcov) return(temp) else return(temp[[1]])},
      predlist,selist,cilblist,ciublist,vcovlist,SIMPLIFY=FALSE)
    # IF PREDICTION ONLY FOR 1 VALUE, ELIMINATE THE FIRST HIERARCHY
    names(fit) <- rownames(X)
    if(length(fit)==1) fit <- fit[[1]]  
  } 
  # SIMPLIFY IF MATRIX AND ONLY 1 PREDICTION
  if(format=="matrix"  && length(predlist)==1L) {
    if(is.matrix(fit)) fit <- fit[1,]
    if(is.list(fit)) fit <- lapply(fit,function(x) x[1,])
  }
#
  fit
}
