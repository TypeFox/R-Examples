summary.mitml <- function(object, n.Rhat=3, goodness.of.appr=FALSE,...){
# summary method for objects of class "mitml"

  # .SDprop <- mitml:::.SDprop
  # .GelmanRubin <- mitml:::.GelmanRubin

  inc <- object$data
  ngr <- length(unique(attr(object$data,"group")))
  prm <- object$par.imputation

  # percent missing
  mdr <- sapply(inc, FUN=function(x){mean(is.na(x))})
  mdr[] <- sprintf(mdr*100,fmt="%.1f")
  mdr <- gsub("^0.0$","0",mdr)

  # convergence for imputation phase
  conv <- NULL
  iter <- dim(prm[[1]])[3]
  Rhat <- ifelse(is.null(n.Rhat), FALSE, n.Rhat >= 2)

  SDprop <- goodness.of.appr
  if(Rhat|SDprop){

    conv <- list(beta=NULL,psi=NULL,sigma=NULL)
    for(pp in c("beta","psi","sigma")){

      ni <- dim(prm[[pp]])[1]
      nj <- dim(prm[[pp]])[2]
      nl <- dim(prm[[pp]])[4]
      cmat <- matrix(NA_real_, ni*nj*nl, 3+Rhat+SDprop)
      cmat[,1] <- rep(1:ni,nj*nl)
      cmat[,2] <- rep(1:nj,each=ni,times=nl)
      cmat[,3] <- rep(1:nl,each=ni*nj)
      colnames(cmat) <- c("i1","i2","grp",if(Rhat) "Rhat",if(SDprop) "SDprop")

      for(ll in 1:nl){ # by group

        lind <- cmat[,3]==ll
        chains <- matrix(prm[[pp]][,,,ll], ni*nj, iter)
        # potential scale reduction (Rhat)
        if(Rhat) cmat[lind,"Rhat"] <- .GelmanRubin(chains,n.Rhat)
        # goodness of approximation
        if(SDprop) cmat[lind,"SDprop"] <- .SDprop(chains)
      
      }
      conv[[pp]] <- cmat
    }

  attr(conv,"stats") <- c("Rhat","SDprop")[c(Rhat,SDprop)]
  }

  smr <- list(
    call=object$call,
    model=object$model,
    prior=object$prior,
    iter=object$iter,
    ngr=ngr,
    missing.rates=mdr,
    conv=conv
  )
  
  class(smr) <- "mitml.summary"
  smr
}

