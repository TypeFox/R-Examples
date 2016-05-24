getsm <-
function(object, t.pval=0.05, wald.pval=t.pval,
  vcov.type=NULL, do.pet=TRUE,
  ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL,pval=0.025),
  normality.JarqueB=NULL, info.method=c("sc", "aic", "hq"),
  keep=NULL, include.gum=FALSE, include.empty=FALSE,
  max.regs=NULL, zero.adj=NULL, vc.adj=NULL, verbose=TRUE,
  print.searchinfo=TRUE, estimate.specific=TRUE, plot=TRUE,
  alarm=FALSE)
{
  ### ARGUMENTS: ###########

  info.method <- match.arg(info.method)

  if(is.null(vcov.type)){
    vcov.type <- object$aux$vcov.type
  }else{
    types <- c("ordinary", "white", "newey-west")
    which.type <- charmatch(vcov.type, types)
    vcov.type <- types[which.type]
  }
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
  if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }
  tol <- object$aux$tol
  LAPACK <- object$aux$LAPACK

  ### INITIALISE ##########

  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  notes <- list()
  spec <- list()
  spec.results <- NULL

##in the future, check object$aux$vX instead?:

  ##check if variance equation:
  if(is.null(object$variance.results)){
    var.spec.chk <- FALSE
  }else{var.spec.chk <- TRUE}

  ## GUM: ##################

  if(is.null(object$aux$mX)){ stop("Mean equation empty") }

  ##deletable, non-deletable regressors, re-organise:
  keep.n <- length(keep)
  gum <- 1:object$aux$mXncol
  delete <- setdiff(gum, keep)
  delete.n <- length(delete)
  if(delete.n > 0){mXdel <- cbind(object$aux$mX[,delete])}else{mXdel <- NULL}
  if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(object$aux$mX[,keep])}
  mXadj <- cbind(mXdel,mXndel)

  ##estimate GUM:
  est <- ols(object$aux$y, mXadj, tol = object$aux$tol,
    LAPACK=LAPACK, method=2)
  fit <- as.vector(mXadj%*%cbind(est$coefficients))
  resids2 <- coredata(object$resids)^2
  d.f. <- object$aux$y.n - object$aux$mXncol
  sumResids2 <- sum(resids2)
  sigma2 <- sumResids2/d.f.

  ##estimate s.e.; compute t-stats. and p-vals.:
  if(vcov.type == "ordinary"){
    varcovmat <- sigma2*est$xtxinv
    coef.var <-as.vector(diag(varcovmat))
    s.e. <- sqrt(coef.var)
  }
  if(vcov.type == "white"){
      omega.hat <- crossprod(mXadj, mXadj*resids2)
#OLD:
#    matResids2 <- matrix(0, object$aux$y.n, object$aux$y.n)
#    diag(matResids2) <- resids2
#    omega.hat <- t(mXadj)%*%matResids2%*%mXadj
    varcovmat <- est$xtxinv %*% omega.hat %*% est$xtxinv
    coef.var <- as.vector(diag(varcovmat))
    s.e. <- sqrt(coef.var)
  }
  if(vcov.type == "newey-west"){
    iL <- round(object$aux$y.n^(1/4), digits=0)
    vW <- 1 - 1:iL/(iL+1)
    vWsqrt <- sqrt(vW)
    mXadjj <- object$resids*mXadj
    mS0 <- crossprod(mXadjj)

    mSum <- 0
    for(l in 1:iL){
      mXadjjw <- mXadjj*vWsqrt[l]
      mXadjjwNo1 <- mXadjjw[-c(1:l),]
      mXadjjwNo2 <- mXadjjw[-c(c(object$aux$y.n-l+1):object$aux$y.n),]
      mSum <- mSum + crossprod(mXadjjwNo1, mXadjjwNo2) + crossprod(mXadjjwNo2, mXadjjwNo1)
    }

    omega.hat <- mS0 + mSum
    varcovmat <- est$xtxinv %*% omega.hat %*% est$xtxinv
    coef.var <- as.vector(diag(varcovmat))
    s.e. <- sqrt(coef.var)
  }
  t.stat <- est$coefficients/s.e.
  p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

  ##diagnostics:
  diagnostics.table <- object$diagnostics[1:3,]

  ##Ljung-Box test of zhat:
  if(is.null(ar.LjungB)){ ar.gum.chk <- 1 }else{

    ##check if new test is needed:
    if(ar.LjungB[1]==object$aux$qstat.options[1]){
      ar.LjungB.new <- FALSE
    }else{ ar.LjungB.new <- TRUE }

    ##if new test needed:
    if(ar.LjungB.new){
      ar.LjungBox <- Box.test(coredata(na.trim(object$resids.std)),
        lag = ar.LjungB[1], type="L")
      if(ar.LjungBox$p.value <= ar.LjungB[2]){
        ar.gum.chk <- 0
      }else{ ar.gum.chk <- 1 }
      rownames(diagnostics.table)[1] <- paste("Ljung-Box AR(",
        ar.LjungB[1], ")", sep="")
      diagnostics.table[1,1] <- ar.LjungBox$statistic
      diagnostics.table[1,2] <- ar.LjungB[1]
      diagnostics.table[1,3] <- ar.LjungBox$p.value
    }

    ##new test not needed:
    if(!ar.LjungB.new){
      if(diagnostics.table[1,3] <= ar.LjungB[2]){
        ar.gum.chk <- 0
      }else{ ar.gum.chk <- 1 }
    }

  } #end if(is.null(..))

  #Ljung-Box test of zhat^2:
  if(is.null(arch.LjungB)){ arch.gum.chk <- 1 }else{

    ##check if new test is needed:
    if(arch.LjungB[1]==object$aux$qstat.options[2]){
      arch.LjungB.new <- FALSE
    }else{ arch.LjungB.new <- TRUE }

    ##if new test needed:
    if(arch.LjungB.new){
      arch.LjungBox <- Box.test(coredata(na.trim(object$resids.std))^2,
        lag = arch.LjungB[1], type="L")
      if(arch.LjungBox$p.value <= arch.LjungB[2]){
        arch.gum.chk <- 0
      }else{ arch.gum.chk <- 1 }
      rownames(diagnostics.table)[2] <- paste("Ljung-Box ARCH(",
        arch.LjungB[1], ")", sep="")
      diagnostics.table[2,1] <- arch.LjungBox$statistic
      diagnostics.table[2,2] <- arch.LjungB[1]
      diagnostics.table[2,3] <- arch.LjungBox$p.value
    }

    ##new test not needed:
    if(!arch.LjungB.new){
      if(diagnostics.table[2,3] <= arch.LjungB[2]){
        arch.gum.chk <- 0
      }else{ arch.gum.chk <- 1 }
    }

  } #end if(is.null(arch..

  ## Jarque-Bera test of zhat:
  if(is.null(normality.JarqueB)){ normality.gum.chk <- 1 }else{
    if(diagnostics.table[3,3] <= normality.JarqueB[1]){
      normality.gum.chk <- 0
    }else{ normality.gum.chk <- 1 }
  }#end if(is.null(normality..

  ## results:
  if(verbose){
    tmp <- rep(0,object$aux$mXncol)
    if(!is.null(keep)){ tmp[keep] <- 1 }
    tmpdf <- cbind(tmp, object$mean.results)
    tmp <- 1:object$aux$mXncol
    tmpdf <- cbind(tmp, tmpdf)
    colnames(tmpdf)[1:2] <- c("reg.no", "keep")
    out$gum.mean <- tmpdf
    out$gum.variance <- object$variance.results
    out$gum.diagnostics <- diagnostics.table
  } #end if(verbose)

  ## if GUM passes diagnostic checks:
  gum.chk <- ar.gum.chk*arch.gum.chk*normality.gum.chk
  if(gum.chk != 0){

    spec[[1]] <- spec.gum <- gum

    #specification results
    info.results <- info.criterion(object$logl, object$aux$loge2.n,
      object$aux$mXncol, method = info.method)
    spec.results <- rbind( c(info.results$value, object$logl,
      info.results$n, info.results$k) )
    col.labels <- c(paste("info(", info.method, ")", sep=""),
      "logl", "n", "k")
    row.labels <- c("spec 1 (gum):")

    #record data for tests against gum:
    gum.regs <- c(delete, keep)
    gum.coefs <- object$mean.results[gum.regs,1]
    gum.varcovmat <- varcovmat

  }else{
    notes <- c(notes, c("MGUM does not pass one or more diagnostic checks"))
  }

  ## FUTURE: ADD 1-CUT MODEL #########

  ## EMPTY MODEL: ################

  if( gum.chk !=0 && delete.n>0 && include.empty==TRUE ){

    ## DO NOT do pet in order to enable reality check!

    ## estimate model:
    if(!is.null(mXndel)){
      est <- ols(object$aux$y, mXndel, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, method=2)
      fit <- as.vector(mXndel%*%cbind(est$coefficients))
    }else{
      fit <- rep(0, object$aux$y.n)
    }
    resids <- as.vector(object$aux$y) - fit
    resids2 <- resids^2
    mXndel.k <- NCOL(mXndel)
    d.f. <- object$aux$y.n - mXndel.k
    sumResids2 <- sum(resids2)
    sigma2 <- sumResids2/d.f.

    #make resids.adj, sigma2.fit, zhat:
    if(var.spec.chk){
      residsAdj <- zoo(resids, order.by=object$aux$y.index)
      est.var <- arx(residsAdj, vc=object$aux$vc,
#OLD:
#      est.var <- arx(resids, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
        zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
        tol=object$aux$tol, LAPACK=object$aux$LAPACK,
        verbose=TRUE, plot=FALSE)
      sigma2.fit <- coredata(na.trim(est.var$var.fit))
      resids.adj <- resids[c(object$aux$y.n-object$aux$loge2.n+1):object$aux$y.n]
      resids.adj.n <- length(resids.adj)
#OLD (ERROR):
#      resids.adj.n <- object$aux$y.n
      zhat <- coredata(na.trim(est.var$resids.std))
    }else{
      resids.adj <- resids
      resids.adj.n <- length(resids.adj)
#OLD (ERROR):
#      resids.adj.n <- object$aux$y.n
      sigma2.fit <- rep(sigma2, resids.adj.n)
      zhat <- resids/sqrt(sigma2)
    } #end if(var.spec.chk)

    diagnostics.chk <- diagnostics(zhat, s2=1,
      ar.LjungB=ar.LjungB, arch.LjungB=arch.LjungB,
      normality.JarqueB=normality.JarqueB, verbose=FALSE)

    ##if empty model passes diagnostic checks:
    if(diagnostics.chk){

      ##add empty to spec:
      spec[[length(spec)+1]] <- if(is.null(keep)){0}else{keep}

      ##specification results
      logl <- -resids.adj.n*log(2*pi)/2 - sum(log(sigma2.fit))/2 - sum(zhat^2)/2
      info.results <- info.criterion(logl, resids.adj.n, keep.n,
        method = info.method)

      ##add empty to spec:
      spec.results <- rbind(spec.results,
        c(info.results$value, logl, info.results$n,
        info.results$k))
      row.labels <- c(row.labels,
        paste("spec ", length(spec), " (empty):", sep=""))

    }else{
        notes <- c(notes, c("Empty mean model does not pass one or more diagnostic checks"))
    } #end if(empty passes diagnostics==TRUE){..}else{..}

  } #end if(include empty model==TRUE)

## MULTI-PATH SEARCH: #################

#future:
#regs <- list()
#regs.info <- list()

#the following should probably be moved to just before
#"single path search", that is, just beneath "prepare
#single-path search":
#regs.current.path <- list() #the regressions of the current path
#regs.info.current.path <- list() #the regression info of current path
#regression info: list(which.path=??, where.in.path=??)

insig.regs <- NULL
paths <- list()
if( gum.chk!=0 && delete.n>1 ){

  ## paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  n.paths <- length(insig.regs)

  ## if paths = 0:
  if(n.paths == 0){
    notes <- c(notes, c("All regressors significant in GUM mean equation"))
  }

  ## if paths > 0:
  if(n.paths > 0){

    ## paths:
    for(i in 1:n.paths){

      ## print searchinfo:
      if(print.searchinfo){
        cat("Searching path no. ", i, " out of ",
          n.paths, "\n", sep="")
      }

      ## prepare single-path search:
      path <- insig.regs[i]
      delete.adj <- setdiff(delete, insig.regs[i])
      keep.adj <- as.numeric(keep)

      ## single-path search of path i:
      for(j in 1:max.regs){

        ## matrices:
        mXdell <- if(length(delete.adj)==0){NULL}else{object$aux$mX[,delete.adj]}
        mXndell <- if(is.null(keep.adj)){NULL}else{object$aux$mX[,keep.adj]}

        ## estimate model:
        mXadj <- cbind(mXdell,mXndell)
        if(!is.null(mXadj)){
          est <- ols(object$aux$y, mXadj, tol=object$aux$tol,
            LAPACK=object$aux$LAPACK, method=2)
          fit <- as.vector(mXadj%*%cbind(est$coefficients))
        }else{
          fit <- rep(0, object$aux$y.n)
        }
        resids <- as.vector(object$aux$y) - fit
        resids2 <- resids^2
        mXadj.k <- NCOL(mXadj)
        d.f. <- object$aux$y.n - mXadj.k
        sumResids2 <- sum(resids2)
        sigma2 <- sumResids2/d.f.

        ## make resids.adj, sigma2.fit, zhat:
        if(var.spec.chk){
          residsAdj <- zoo(resids, order.by=object$aux$y.index)
#OLD:
#          est.var <- arx(resids, vc=object$aux$vc,
          est.var <- arx(residsAdj, vc=object$aux$vc,
            arch=object$aux$arch, asym=object$aux$asym,
            log.ewma=object$aux$log.ewma,
            vxreg=object$aux$vxreg,
            zero.adj=object$aux$zero.adj,
            vc.adj=object$aux$vc.adj, tol=object$aux$tol,
            LAPACK=object$aux$LAPACK, verbose=TRUE, plot=FALSE)
          sigma2.fit <- coredata(na.trim(est.var$var.fit))
          resids.adj <- resids[c(object$aux$y.n-object$aux$loge2.n+1):object$aux$y.n]
          resids.adj.n <- length(resids.adj)
#OLD (ERROR):
#          resids.adj.n <- object$aux$y.n
          zhat <- coredata(na.trim(est.var$resids.std))
        }else{
          resids.adj <- resids
          resids.adj.n <- length(resids.adj)
#OLD (ERROR):
#          resids.adj.n <- object$aux$y.n
          sigma2.fit <- rep(sigma2, resids.adj.n)
          zhat <- resids/sqrt(sigma2)
        } #end if(var.spec.chk)

        diagnostics.chk <- diagnostics(zhat, s2=1,
          ar.LjungB=ar.LjungB, arch.LjungB=arch.LjungB,
          normality.JarqueB=normality.JarqueB, verbose=FALSE)

        ## if diagnostics.chk fails (i.e. FALSE),
        ## then move path[length(path)] to keep.adj
        if(!diagnostics.chk){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*c(-1))
          next #next j
        }

        #if diagnostics are ok (i.e. TRUE):
        if(diagnostics.chk){

          ## stop if no more deletable regressors:
          if(length(delete.adj)==0){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          if(!is.null(mXadj)){

            ## estimate s.e.; compute t-stats. and p-vals.:
            if(vcov.type == "ordinary"){
              coef.var <-as.vector(sigma2*diag(est$xtxinv))
              s.e. <- sqrt(coef.var)
            } #end "ordinary"

            if(vcov.type == "white"){
              omega.hat <- crossprod(mXadj, mXadj*resids2)
#OLD:
#              matResids2 <- matrix(0, object$aux$y.n, object$aux$y.n)
#              diag(matResids2) <- resids2
#              omega.hat <- t(mXadj)%*%matResids2%*%mXadj
              varcovmat <- est$xtxinv %*% omega.hat %*% est$xtxinv
              coef.var <- as.vector(diag(varcovmat))
              s.e. <- sqrt(coef.var)
            } #end "white"

            if(vcov.type == "newey-west"){
              iL <- round(object$aux$y.n^(1/4), digits=0)
              vW <- 1 - 1:iL/(iL+1)
              vWsqrt <- sqrt(vW)
              mXadjj <- resids*mXadj
              mS0 <- crossprod(mXadjj)

              mSum <- 0
              for(l in 1:iL){
                mXadjjw <- mXadjj*vWsqrt[l]
                mXadjjwNo1 <- mXadjjw[-c(1:l),]
                mXadjjwNo2 <- mXadjjw[-c(c(object$aux$y.n-l+1):object$aux$y.n),]
                mSum <- mSum + crossprod(mXadjjwNo1, mXadjjwNo2) + crossprod(mXadjjwNo2, mXadjjwNo1)
              }

              omega.hat <- mS0 + mSum
              varcovmat <- est$xtxinv %*% omega.hat %*% est$xtxinv
              coef.var <- as.vector(diag(varcovmat))
              s.e. <- sqrt(coef.var)
            } #end newey-west

            ## t-tests:
            t.stat <- est$coefficients/s.e.
            p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

          } #end if(!is.null(mXadj))

          ## if any p-value > t.pval:
          if(sum(p.val[1:c(length(delete.adj))] > t.pval) > 0){

            reg.no <- which.max(p.val[1:I(length(delete.adj))])

            ## do pet test:
            if(do.pet){
              deleted <- setdiff(delete, delete.adj[-reg.no])
              n.deleted <- length(deleted)
              mR <- NULL #initiate restriction matrix
              for(k in 1:object$aux$mXncol){
                if(gum.regs[k] %in% deleted){
                  mR <- rbind(mR, c(rep(0,I(k-1)), 1, rep(0, I(object$aux$mXncol-k) )))
                } #end if(gum.regs[k}..)
              } #end for(k in ..)
              mRestq <- mR%*%cbind(gum.coefs)
              wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=object$aux$tol)%*%mRestq
              pet.chk <- as.logical(wald.pval < pchisq(wald.stat, n.deleted, lower.tail = FALSE))
            }else{
              pet.chk <- TRUE
            } #end if(do.pet)else..

            ## delete regressor if(pet.chk), else move to keep:
            if(pet.chk){
              path <- union(path, delete.adj[reg.no])
              delete.adj <- delete.adj[-reg.no]
            }else{
              path <- union(path, delete.adj[reg.no]*I(-1))
              keep.adj <- union(delete.adj[reg.no], keep.adj)
              delete.adj <- delete.adj[-reg.no]
            } #end if(pet.chk)else{..}

          }else{
            spec.adj <- union(delete.adj, keep.adj)
            break
          } #end if..else.. any p-value > t.pval

        } #end if diagnostics are ok

      } #### end single-path search: for(j in..

      #it is probably at this point that I should introduce
      #the 'bookkeeping' with respect to regs, regs.info,
      #regs.current.path and regs.info.current.path

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj is already in spec:
      if(length(spec.adj)==0){spec.adj <- 0} #check if completely empty
      for(l in 1:length(spec)){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break} #stop for(l in..)
      }

      #if spec.adj not in spec (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to spec:
        spec[[length(spec)+1]] <- spec.adj

        #specification results
        if(spec.adj[1]==0){
          n.spec.adj <- 0
        }else{
          n.spec.adj <- length(spec.adj)
        }
        logl <- -resids.adj.n*log(2*pi)/2 - sum(log(sigma2.fit))/2 - sum(zhat^2)/2
        info.results <- info.criterion(logl, resids.adj.n,
          n.spec.adj, method=info.method)

        #add terminal to spec.results:
        spec.results <- rbind(spec.results,
          c(info.results$value, logl, info.results$n,
          info.results$k))
        row.labels <- c(row.labels, paste("spec ", length(spec), ":", sep=""))

      } #end if(chk.spec==FALSE)

    } #end multi-path search: for(i in 1:n.paths) loop

  } #end if paths > 0
} #end if( gum.chk!=0 && delete.n>1 )

  ## FIND THE BEST MODEL: ########################

  #future?: check first if spec results is empty, then add
  #gum to it if it is?

  if(!is.null(spec.results)){

    J <- 1:NROW(spec.results)
    models <- cbind(J, spec.results)
    colnames(models) <- NULL

    #find best model and check for several minimums:
    if(include.gum){
      min.value <- min(models[,2])
      where <- which(min.value==models[,2])
    }else{
      if(length(spec)==1){
        where <- 1
#        notes <- c(notes, "MGUM added to terminals since there were no other terminals")
      }else{
        min.value <- min(models[-1,2])
        where <- which(min.value==models[-1,2]) + 1
      } #end if(length(spec)==1)
    } #end if(include.gum)..
    if(length(where)>1){ notes <- c(notes, "Several terminal specifications attain the minimum information criterion") }
    best.spec <- spec[[where[1]]] #winner

  } #end if(!is.null(spec.results))

  ## OUTPUT ################################

  out$keep <- keep
  out$insigs.in.gum <- insig.regs

  ##if no search has been undertaken:
  if(is.null(spec.results)){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }

  ##if search has been undertaken:
  if(!is.null(spec.results)){

    ##terminals results:
    if(length(paths)==0){
      out$paths <- NULL
    }else{ out$paths <- paths }
    out$terminals <- spec
    colnames(spec.results) <- col.labels
    where.empty <- which(spec.results[,"k"]==0)
    if(include.empty==FALSE && length(where.empty) > 0){
      row.labels[where.empty] <- paste("spec ", where.empty,
        " (empty):", sep="")
    }
    rownames(spec.results) <- row.labels
    out$terminals.results <- spec.results

    if(!estimate.specific){
      if(best.spec==0 || is.na(best.spec) || length(best.spec)==0 ){
        out$specific.spec <- NULL
      }else{
        specific <- sort(best.spec)
        names(specific) <- object$aux$mXnames[specific]
        out$specific.spec <- specific
      }
    } #end if(!estimate.specific)

    if(estimate.specific){

#      if(best.spec[1]==0){
#        resids <- zoo(cbind(object$aux$y), order.by=object$aux$y.index)
#        colnames(resids) <- object$aux$y.name
#        specific.mean <- "empty"
#        Rsquared <- 0
#        vcov.mean <- NULL
#      }else{

      ##prepare for estimation:
      yadj <- zoo(cbind(object$aux$y),
        order.by=object$aux$y.index)
      colnames(yadj) <- object$aux$y.name
      specific <- sort(best.spec)
      if(specific[1]==0){
        mXadj <- NULL
      }else{
        mXadj <- cbind(object$aux$mX[,specific])
        colnames(mXadj) <- object$aux$mXnames[specific]
        mXadj <- zoo(mXadj, order.by=object$aux$y.index)
      }
      if(is.null(object$aux$vxreg)){
        vxregAdj <- NULL
      }else{
        vxregAdj <- zoo(object$aux$vxreg,
          order.by=object$aux$y.index)
      }
      if(is.null(ar.LjungB)){
        ar.LjungB <- object$aux$qstat.options[1]
      }
      if(is.null(arch.LjungB)){
        arch.LjungB <- object$aux$qstat.options[2]
      }

      ##estimate specific model:
      est <- arx(yadj, mxreg=mXadj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=vxregAdj,
        zero.adj=object$aux$zero.adj,
        vc.adj=object$aux$vc.adj, vcov.type=vcov.type,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        tol=object$aux$tol, LAPACK=object$aux$LAPACK,
        verbose=TRUE, plot=FALSE)

      ##rename various stuff:
      est$call <- est$date <- NULL
#        where.mean.results <- which(names(est)=="mean.results")
#        if(length(where.mean.results)>0){
#          names(est)[where.mean.results] <- "specific.mean"
#        }
#        where.variance.results <- which(names(est)=="variance.results")
#        if(length(where.variance.results)>0){
#          names(est)[where.variance.results] <- "specific.variance"
#        }
      where.diagnostics <- which(names(est)=="diagnostics")
      est$diagnostics <- est$diagnostics[1:3,]
      if(length(where.diagnostics)>0){
        names(est)[where.diagnostics] <- "specific.diagnostics"
      }
      est$aux$y.name <- object$aux$y.name
      est <- unclass(est)
      names(specific) <- colnames(mXadj)
      out$specific.spec <- specific
      out <- c(out,est)

    } #end if(estimate.specific)

  } #end if(!is.null(spec.results))

  if(length(notes) > 0){ out$notes <- notes }
  out$aux$mXnames.gum <- object$aux$mXnames
  if(is.null(out$aux$vcov.type)){ out$aux$vcov.type <- vcov.type }
  #if(is.null(out$aux$y.n)){ out$aux$y.n <- object$aux$y.n }
  out <- c(list(date=date(), gets.type="getsm"), out)
  out$time.finished <- date()
  class(out) <- "gets"

  if(alarm){ alarm() }
  if(plot){ plot.gets(out) }
  return(out)
}
