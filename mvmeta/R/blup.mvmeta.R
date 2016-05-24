###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
blup.mvmeta <-
function(object, se=FALSE, pi=FALSE, vcov=FALSE, pi.level=0.95,
  format=c("matrix","list"), aggregate=c("stat","y"), ...) {
#
################################################################################
# CREATE DESIGN MATRIX X, S, FITTED VALUES AND RESIDUALS
#
  mf <- model.frame(object)
  X <- model.matrix(object)
  y <- as.matrix(model.response(mf,"numeric"))
  S <- object$S
  offset <- object$offset
  # IN CASE, RE-INSERT MISSING VALUES FROM FITTING (ORIGINAL OBJECT)
  if(!is.null(object$na.action)) {
    X <- napredict(object$na.action,X)
    y <- napredict(object$na.action,y)
    S <- napredict(object$na.action,S)
    offset <- napredict(object$na.action,offset)
  }
  # DEFINE MISSING (THOSE STUDIES WITH ALL MISSING PREDICTED TO NA LATER)
  nay <- is.na(y)
  # SET MISSING RESPONSES TO 0 (TO GET PREDICTION FOR PARTIALLY MISSING)
  y[nay] <- 0
#
################################################################################
# RE-CREATE OBJECTS
#
  format <- match.arg(format,c("matrix","list"))
  aggregate <- match.arg(aggregate,c("stat","y"))
#
  # TRANSFORM X AND y IN LISTS
  Xlist <- lapply(seq(nrow(X)),
    function(i) diag(1,object$dim$k)%x%X[i,,drop=FALSE])
  ylist <- lapply(seq(nrow(X)),function(i) y[i,])
  # PREDICTED VALUES (INCLUDING OFFSET)
  predlist <- lapply(seq(nrow(X)),function(i) {
    pred <- Xlist[[i]]%*%as.vector(object$coefficients)
    if(!is.null(offset)) pred <- pred+offset[i]
    return(pred)})
  # RESIDUALS
  reslist <- mapply(function(y,pred) y-pred,ylist,predlist,SIMPLIFY=FALSE)
#
  # COMPUTE Sigma, WITH VARIANCES CORRESPONDING TO MISSING FILLED IN WITH 10^10
  # NB: FOR COMPLETELY MISSING OBS, SET TO NA LATER THROUGH X
  if(dim(S)[2]==ncol(y)) S <- inputcov(sqrt(S),object$control$Scor)
  Psi <- if(!is.null(object$Psi)) object$Psi else diag(0,object$dim$k)
  Sigmalist <- lapply(seq(nrow(X)),function(i) {
    if(all(nay[i,])) return(diag(1,object$dim$k))
    Si <- diag(10^10,object$dim$k)
    Si[!nay[i,],!nay[i,]] <- xpndMat(S[i,])[!nay[i,],!nay[i,]]
    return(Si+Psi)})
  # COMPUTE ITS INVERSE OF CHOLESKY DECOMPOSITION AND RESIDUALS
  Ulist <- lapply(Sigmalist,chol)
  invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))  
#
################################################################################
#
  # COMPUTE BLUP
  bluplist <- mapply(function(pred,invU,res) {
    blup <- drop(pred + Psi%*%(tcrossprod(invU))%*%res)
    names(blup) <- object$lab$k
    return(blup)},predlist,invUlist,reslist,SIMPLIFY=FALSE)
#
  # COMPUTE VCOV (SUM OF COMPONENTS), SE AND PI
  zvalpi <- qnorm((1-pi.level)/2,lower.tail=FALSE)
  Vlist <- mapply(function(X,invU) {
    V <- X%*%tcrossprod(object$vcov,X) +
      Psi - Psi%*%(tcrossprod(invU)%*%Psi)
    dimnames(V) <- list(object$lab$k,object$lab$k)
    return(V)},Xlist,invUlist,SIMPLIFY=FALSE)
  selist <- lapply(Vlist, function(x) sqrt(diag(x)))
  pilblist <- mapply(function(blup,se) blup-zvalpi*se,
      bluplist,selist,SIMPLIFY=FALSE)
  piublist <- mapply(function(blup,se) blup+zvalpi*se,
      bluplist,selist,SIMPLIFY=FALSE)
#
################################################################################ 
# AGGREGATE AND RETURN
#
  # IF VCOV WHEN MORE THAN 1 OUTCOME, SWITCH TO LIST
  if(vcov && object$dim$k>1) format <- "list"
  # IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
  if(format=="matrix" && object$dim$k==1) {
    blup <- drop(do.call("rbind",bluplist))
    names(blup) <- rownames(X)
    if(se) blup <- cbind(blup,se=drop(do.call("rbind",selist)))
    if(pi) blup <- cbind(blup,pi.lb=drop(do.call("rbind",pilblist)),
      pi.ub=drop(do.call("rbind",piublist)))
    if(vcov) blup <- cbind(blup,vcov=drop(do.call("rbind",Vlist)))
  # IF NO OTHER STAT, FORCE THE AGGREGATION ON PREDICTION
  } else if(format=="matrix" && !se && !pi) {
    blup <- do.call("rbind",bluplist)
    rownames(blup) <- rownames(X)
    # WHEN AGGREGATE ON STAT
  } else if(format=="matrix" && aggregate=="stat") {
    blup <- list(blup=do.call("rbind",bluplist))
    rownames(blup$blup) <- rownames(X)
    if(se) {
      blup$se <- do.call("rbind",selist)
      rownames(blup$se) <- rownames(X)
    }
    if(pi) {
      blup$pi.lb <- do.call("rbind",pilblist)
      blup$pi.ub <- do.call("rbind",piublist)
      rownames(blup$pi.lb) <- rownames(blup$pi.ub) <- rownames(X)
    }
  # WHEN AGGREGATE ON OUTCOME
  } else if (format=="matrix" && aggregate=="y") {
    blup <- list()
    for(j in seq(object$dim$k)) {
      templist <- lapply(seq(bluplist),function(i) {
        blup <- bluplist[[i]][j]
        if(se) blup <- cbind(blup,selist[[i]][j])
        if(pi) blup <- cbind(blup,pilblist[[i]][j],
          piublist[[i]][j])
        return(as.matrix(blup))})
      temp <- do.call("rbind",templist)
      dimnames(temp) <- list(rownames(X),
        c("blup","se","pi.lb","pi.ub")[c(TRUE,se,pi,pi)])
      blup[[j]] <- temp
    }
    names(blup) <- object$lab$k
  # ALL THE OTHER COMBINATIONS, PRODUCE A LIST
  } else {
    blup <- mapply(function(bluparg,searg,pilbarg,piubarg,vcovarg) {
      temp <- list(blup=bluparg)
      if(se) temp$se <- searg
      if(pi) temp$pi.lb <- pilbarg
      if(pi) temp$pi.ub <- piubarg
      if(vcov) temp$vcov <- vcovarg
      # RETURN ONE OBJECT IF ONLY PRED 
      if(se||pi||vcov) {
        return(temp)
      } else return(temp[[1]])},
      bluplist,selist,pilblist,piublist,Vlist,SIMPLIFY=FALSE)
    names(blup) <- rownames(X)
  }
#
  blup
}
