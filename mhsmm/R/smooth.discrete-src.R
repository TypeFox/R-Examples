##
## Smooth discrete time series using a hidden Markov model
##

smooth.discrete <- function(y, init=NULL, trans=NULL, parms.emission=.5, method="viterbi",
                            details=0,...){
  
  method <- match.arg(method,c("viterbi","smoothed"))

  ## Initial distribution
  ##
  if (is.null(init)){
    init <- as.numeric(table(y))
    init <- init/sum(init)
  }
  J    <- length(init)
  if (details>0){
    cat(sprintf("n.states   : %i\n", J))
    cat(sprintf("init       : %s\n", toString(round(init,4))))
  }

  ## Transition matrix
  if (is.null(trans)){
    ttt<-table(y[-length(y)],y[-1])
    P  <- sweep(ttt, 1, rowSums(ttt), "/")
  } else {
    if (is.matrix(trans)){
      P <- trans
    } else {
      P    <- createTransition(trans,J)
    }
  }
  if (details>0){
    cat(sprintf("transition :\n"))
    print(P)
  }

  if (is.null(parms.emission)){
    BB    <- createTransition(0.8,J)
  } else {
    if (is.matrix(parms.emission)){
      BB <- parms.emission
    } else {
      BB    <- createTransition(parms.emission,J)
    }
  }
  B    <- list(pmf=BB)
  
  if (details>0){
    cat(sprintf("emission   :\n"))
    print(BB)
  }

  init.spec <- hmmspec(init,trans=P,parms.emission=B,dens.emission=.dpmf)
  train   <- list(s=y, x=y, N=length(y))
  
  hmm.obj <- hmmfit(train, init.spec,mstep=.mstep.pmf,...)
  pred <- predict(hmm.obj, train, method=method)
  ans <- list(s=pred$s, model=hmm.obj, data=train, initial=list(init=init, trans=P, parms.emission=B))
  class(ans) <- "smoothDiscrete"
  ans
}

summary.smoothDiscrete <- function(object, ...){
  summary(object$model)
}

print.smoothDiscrete <- function(x,...){
  cat("A 'smoothDiscrete' object\n")
  print(str(x,1))
  return(invisible(x))
}

predict.smoothDiscrete <- function(object, x, method="viterbi", ...){

  method <- match.arg(method,c("viterbi","smoothed"))

  if (missing(x)){
    predict(object$model, x=object$data,method=method)
  } else {
    xxx   <- list(s=x, x=x, N=length(x))
    predict(object$model, xxx)
  }
}

createTransition <- function(Pvec, J){
  Pvec <- rep(Pvec, J)[1:J]
  
  PP <- lapply(Pvec, function(xx)
               rep((1-xx)/(J-1),J)
               )
  PP <- do.call(rbind,PP)
  diag(PP) <- Pvec
  PP
}


##
## Auxillary functions for smooth.discrete
##

.dpmf <- function(x,j,model) {
   ret <- model$parms.emission$pmf[j,x]
   ret[is.na(ret)]=1
   ret
  }

.mstep.pmf <- function(x,wt) {
   ans <- matrix(ncol=ncol(wt),nrow=ncol(wt))
   for(i in 1:ncol(wt))
     for(j in 1:ncol(wt))
       ans[i,j] <- sum(wt[which(x[!is.na(x)]==j),i])/sum(wt[!is.na(x),i])
   list(pmf=ans)
}


