refit <-
  ## Short form for generic function for refitting model
  function(object, newpatterns,newfreq, useinit=FALSE,  ...) UseMethod("refit")

`refit.randomLCA` <-
  function(object,newpatterns,newfreq=NULL,useinit=FALSE,...) {
    if (!inherits(object, "randomLCA"))
      stop("Use only with 'randomLCA' objects.\n")
    if (is.null(newfreq)) {
      pats <- apply(as.matrix(newpatterns), 1, function(x) {paste(ifelse(is.na(x),"N",x),collapse="")})
      tpats <- table(pats)
      newfreq <- as.numeric(tpats)
      new2patterns <- unlist(strsplit(names(tpats),split=""))
      new2patterns <- ifelse(new2patterns=="N",NA_character_,new2patterns)
      new2patterns <- as.data.frame(matrix(as.numeric(new2patterns),byrow=TRUE,ncol=dim(newpatterns)[2]))
      if (is.null(names(newpatterns))) names(new2patterns) <- paste("X",1:dim(newpatterns)[2],sep="")
      else names(new2patterns) <- names(newpatterns)
      newpatterns <- new2patterns
    }
    else {
      # check that newfreq doesn't contain missing
      if (any(is.na(newfreq))) stop("newfreq cannot contain missing values")
      # remove any observations with newfrequency of zero
      newpatterns <- newpatterns[newfreq!=0,]
      newfreq <- newfreq[newfreq!=0]
    }
    if (useinit) {
      if (!object$random) newfit <- fitFixed(newpatterns,newfreq,object$outcomep,object$classp,
                                             object$nclass,calcSE=FALSE,justEM=FALSE,object$probit,object$penalty,verbose=FALSE)
      else {
        if (!object$level2) newfit <- fitAdaptRandom(newpatterns,newfreq,
                                                     nclass=object$nclass,calcSE=FALSE,initoutcomep=object$outcomep,
                                                     initclassp=object$classp,
                                                     initlambdacoef=object$lambdacoef,
                                                     blocksize=object$blocksize,
                                                     gh=norm.gauss.hermite(object$quadpoints),
                                                     constload=object$constload,probit=object$probit,byclass=object$byclass,
                                                     qniterations=object$qniterations,penalty=object$penalty,verbose=FALSE)
        else newfit <- fitAdaptRandom2(newpatterns,newfreq,
                                       nclass=object$nclass,calcSE=FALSE,initoutcomep=object$outcomep,
                                       initclassp=object$classp,
                                       initlambdacoef=object$lambdacoef,
                                       initltaucoef=object$ltaucoef,
                                       level2size=object$level2size,
                                       constload=object$constload,
                                       gh=norm.gauss.hermite(object$quadpoints),
                                       probit=object$probit,byclass=object$byclass,
                                       qniterations=object$qniterations,penalty=object$penalty,verbose=FALSE)
      }
      newfit$call <- object$call
      newfit$nclass <- object$nclass
      newfit$random <- object$random
      newfit$constload <- object$constload
      newfit$level2 <- object$level2
      newfit$byclass <- object$byclass
      newfit$probit <- object$probit
      newfit$quadpoints <- object$quadpoints
      newfit$blocksize <- object$blocksize
      newfit$level2size <- object$level2size
      newfit$patterns <- object$patterns
      newfit$notrials <- object$notrials
      newfit$freq <- object$freq
      newfit$qniterations <- object$qniterations
      newfit$penalty <- object$penalty
      class(newfit) <- "randomLCA"
      } else newfit <- randomLCA(newpatterns,freq=newfreq,nclass=object$nclass,calcSE=TRUE,notrials=object$notrials,
             random=object$random,byclass=object$byclas,quadpoints=object$quadpoints,constload=object$constload,
             blocksize=object$blocksize,
             level2=object$level2,probit=object$probit,level2size=object$level2size,
             qniterations=object$qniterations,penalty=object$penalty,verbose=FALSE)
    return(newfit)    
  }

