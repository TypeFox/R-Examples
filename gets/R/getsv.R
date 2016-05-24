getsv <-
function(object, t.pval=0.05, wald.pval=t.pval,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025),
  normality.JarqueB=NULL, info.method=c("sc", "aic", "hq"),
  keep=c(1), include.gum=FALSE, include.empty=FALSE,
  max.regs=NULL, zero.adj=NULL, vc.adj=NULL, verbose=TRUE,
  print.searchinfo=TRUE, estimate.specific=TRUE, plot=TRUE,
  alarm=FALSE)
{
  ### ARGUMENTS ###########

  info.method <- match.arg(info.method)
  vc=TRUE #obligatory
  vcov.type <- "ordinary" #obligatory
  tol <- object$aux$tol
  LAPACK <- object$aux$LAPACK

  ##zoo and NA related:
  e <- object$resids #should not contain NAs
  e.index <- index(e) #use object$aux$y.index instead?
  e <- coredata(e)
  e.n <- length(e) #use object$aux$y.n instead?
  eadj <- e[c(e.n-object$aux$loge2.n+1):e.n] #Note: log(eadj^2)=loge2
  eadj.n <- length(eadj)
  eadj.index <- e.index[c(e.n-object$aux$loge2.n+1):e.n]

#Do I need this?:
#  ##arch arguments:
#  arch <- object$aux$arch
#  asym <- object$aux$asym
#  log.ewma <- object$aux$log.ewma
#
#Do I need this?:
#  vxreg <- object$aux$vxreg
#  if(!is.null(vxreg)){
#    vxreg <- as.zoo(cbind(vxreg))
#    vxreg.names <- colnames(vxreg)
#    vxreg <- window(vxreg, start=e.index[1], end=e.index[e.n])
#    vxreg <- cbind(coredata(vxreg))
#    colnames(vxreg) <- vxreg.names
#  }

  ##diagnostics options, max.regs:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
  if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }

  ##zero-handling:
  if(is.null(zero.adj)){ zero.adj <- object$aux$zero.adj }
  if(is.null(vc.adj)){ vc.adj <- object$aux$vc.adj }

  ##tolerancy for non-invertibility of matrices:
  if(is.null(tol)){ tol <- object$aux$tol }
  if(is.null(LAPACK)){ LAPACK <- object$aux$LAPACK }

  ### INITIALISE ##########

  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  notes <- list()
  spec <- list()
  spec.results <- NULL

  ## REGRESSAND AND REGRESSORS: #####

#  mXorig <- regs.var(e, vc=1, arch=arch, log.ewma=log.ewma,
#    asym=asym, vxreg=vxreg, zero.adj=zero.adj)
#  mX.colnames <- colnames(mXorig[,-1])
#  if(!is.null(log.ewma)){ewma.chk <- 1}else{ewma.chk <- 0}
#  t1 <- max(arch, asym, ewma.chk) + 1
#  eadj <- e[t1:e.n]
#  logep <- mXorig[t1:e.n, 1]
#  mX <- cbind(mXorig[t1:e.n, -1])
  loge2 <- object$aux$loge2
  mX <- object$aux$vX

  ## GUM: #########################

  if(!(1 %in% keep)){
    notes <- c(notes, "Regressor 1 set to 'keep'")
    keep <- union(keep, 1)
  }
  keep.n <- length(keep)
  gum.n <- object$aux$vXncol #OLD: NCOL(mX)
  gum <- 1:gum.n
  delete <- setdiff(gum, keep)
  delete.n <- length(delete)

  ## deletable and non-deletable regressors:
  if(delete.n > 0){mXdel <- cbind(mX[,delete])}else{mXdel <- NULL}
  if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(mX[,keep])}

  ## estimate GUM:
  mXadj <- cbind(mXdel,mXndel)
  est <- ols(loge2, mXadj, tol = tol, LAPACK=LAPACK,
    method=2)
  fit <- as.vector(mXadj%*%cbind(est$coefficients))
  ustar <- as.vector(loge2) - fit
  ustar2 <- ustar^2
  mXadj.n <- NROW(mXadj)
  mXadj.k <- NCOL(mXadj)
  d.f. <- mXadj.n - mXadj.k
  sumUstar2 <- sum(ustar2)
  sigma2 <- sumUstar2/d.f.

  ## estimate s.e.; compute t-stats. and p-vals.:
  if(vcov.type == "ordinary"){
    varcovmat <- sigma2*est$xtxinv
    coef.var <-as.vector(diag(varcovmat))
    s.e. <- sqrt(coef.var)
  }
  t.stat <- est$coefficients/s.e.
  p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

  ## make standardised residuals {z_hat}:
  Elnz2 <- -log(mean(exp(ustar)))
  vconstadj <- -Elnz2
  sigma2adj <- exp(fit + vconstadj)

  zhat <- eadj/sqrt(sigma2adj)
  out$gum.resids.std <- zoo(c(rep(NA, c(e.n-length(zhat))),
    zhat), order.by=e.index)

  ## adjust variance constant and its p-value:
  if(vc.adj==TRUE){
    where.const <- which(union(delete,keep)==1)
    waldstat <- ((est$coefficients[where.const] + vconstadj)^2)/s.e.[where.const]^2
    p.val[where.const] <- pchisq(waldstat, 1, lower.tail = FALSE)
    est$coefficients[where.const] <- est$coefficients[where.const] + vconstadj
  }

  ## diagnostics:
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

  #Jarque-Bera test of zhat:
  if(is.null(normality.JarqueB)){ normality.gum.chk <- 1 }else{
    if(diagnostics.table[3,3] <= normality.JarqueB[1]){
      normality.gum.chk <- 0
    }else{ normality.gum.chk <- 1 }
  }#end if(is.null(normality..


#OLD:
#  #Ljung-Box test for serial correlation in {z_t}:
#  if(!is.null(ar.LjungB)){
#    ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
#    if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.gum.chk <- 0}else{ar.gum.chk <- 1}
#    if(verbose == TRUE){
#      diagnostics.table[1,1] <- ar.LjungBox$statistic
#      diagnostics.table[1,2] <- ar.LjungBox$p.value}
#  }else{ar.gum.chk <- 1}
#
#  #Ljung-Box test for arch in {z_t^2}
#  if(!is.null(arch.LjungB)){
#    arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
#    if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.gum.chk <- 0}else{arch.gum.chk <- 1}
#    if(verbose == TRUE){
#      diagnostics.table[2,1] <- arch.LjungBox$statistic
#      diagnostics.table[2,2] <- arch.LjungBox$p.value}
#  }else{arch.gum.chk <- 1}

  #results:
  if(verbose){
    tmp <- rep(0,object$aux$vXncol)
    tmp[keep] <- 1
    tmpdf <- cbind(tmp, object$variance.results)
    tmp <- 1:object$aux$vXncol
    tmpdf <- cbind(tmp, tmpdf)
    colnames(tmpdf)[1:2] <- c("reg.no", "keep")
    out$gum.variance <- tmpdf
    out$gum.diagnostics <- diagnostics.table
#OLD:
#    keep.labels <- c(rep(0,delete.n), rep(1,keep.n))
#    results <- as.data.frame(cbind(c(delete,keep), keep.labels,
#      est$coefficients, s.e., t.stat, p.val))
#    colnames(results) <- c("reg.no", "keep", "coef", "std.error", "t-stat", "p-value")
#    rownames(results) <- mX.colnames[c(delete,keep)] #NULL
#    out$gum.variance <- results
#    out$gum.diagnostics <- diagnostics.table
  } #end if(verbose)

  #if GUM passes diagnostic checks:
  gum.chk <- ar.gum.chk*arch.gum.chk*normality.gum.chk
  if(gum.chk != 0){

    spec[[1]] <- spec.gum <- gum

    #specification results
    logl <- -object$aux$loge2.n*log(2*pi)/2 - sum(log(sigma2adj))/2 - sum(eadj^2/sigma2adj)/2
    info.results <- info.criterion(logl, object$aux$loge2.n,
      object$aux$vXncol, method = info.method)
    spec.results <- rbind( c(info.results$value, logl,
      info.results$n, info.results$k) )
    col.labels <- c(paste("info(", info.method, ")", sep=""),
      "logl", "n", "k")
    row.labels <- c("spec 1 (gum):")

    #record data for Wald-tests against gum:
    gum.regs <- c(delete, keep)
    gum.coefs <- object$variance.results[gum.regs,1]
#    gum.coefs <- est$coefficients
    gum.varcovmat <- varcovmat

  }else{
    notes <- c(notes, c("GUM does not pass one or more diagnostic checks"))
  }

  ##adjust variance constant and its p-value:
  #if(vc.adj==TRUE){
  #  where.const <- which(union(delete,keep)==1)
  #  waldstat <- ((est$coefficients[where.const] + vconstadj)^2)/s.e.[where.const]^2
  #  p.val[where.const] <- pchisq(waldstat, 1, lower.tail = FALSE)
  #  est$coefficients[where.const] <- est$coefficients[where.const] + vconstadj
  #}

  ## FUTURE: ADD 1-CUT MODEL #########

  ## EMPTY MODEL: ################

  if( gum.chk !=0 && delete.n>0 && include.empty==TRUE ){

    #DO NOT do pet in order to enable reality check!

    #estimate model:
    est <- ols(loge2, mXndel, tol = tol, LAPACK=LAPACK,
      method=2)
    fit <- as.vector(mXndel%*%cbind(est$coefficients))
    ustar <- as.vector(loge2) - fit
    ustar.n <- length(ustar)
    mXndel.k <- NCOL(mXndel)
    ustar2 <- ustar^2
    sumUstar2 <- sum(ustar2)
    sigma2 <- sumUstar2/(ustar.n-mXndel.k)

    #make standardised residuals {z_hat}:
    Elnz2 <- -log(mean(exp(ustar)))
    vconstadj <- -Elnz2
    sigma2adj <- exp(fit + vconstadj)
    zhat <- eadj/sqrt(sigma2adj)

    diagnostics.chk <- diagnostics(zhat, s2=1,
      ar.LjungB=ar.LjungB, arch.LjungB=arch.LjungB,
      normality.JarqueB=normality.JarqueB, verbose=FALSE)

#OLD:
#    #Ljung-Box test for serial correlation in {z_t}:
#    if(!is.null(ar.LjungB)){
#      ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
#      if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.bench.chk <- 0}else{ar.bench.chk <- 1}
#    }else{ar.bench.chk <- 1}
#
#    #Ljung-Box test for arch in {z_t}
#    if(!is.null(arch.LjungB)){
#      arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
#      if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.bench.chk <- 0}else{arch.bench.chk <- 1}
#    }else{arch.bench.chk <- 1}

    ##if empty model passes diagnostic checks:
    if(diagnostics.chk){

      ##add empty to spec:
      spec[[length(spec)+1]] <- keep

      ##specification results
      logl <- -eadj.n*log(2*pi)/2 - sum(log(sigma2adj))/2 - sum(eadj^2/sigma2adj)/2
      info.results <- info.criterion(logl, eadj.n, mXndel.k,
        method=info.method)

      ##add empty to spec:
      spec.results <- rbind(spec.results,
        c(info.results$value, logl, info.results$n,
        info.results$k))
      row.labels <- c(row.labels,
        paste("spec ", length(spec), " (empty):", sep=""))

    }else{
      notes <- c(notes, c("Empty log-variance model does not pass one or more diagnostic checks"))
    } #end if(diagnostics.chk)
  } #end if(include.empty)

## MULTI-PATH SEARCH: #################

insig.regs <- NULL
paths <- list()
if( gum.chk!=0 && delete.n>1 ){

  ## paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  n.paths <- length(insig.regs)

  ## if paths = 0:
  if(n.paths == 0){
    notes <- c(notes, c("All regressors significant in GUM variance equation"))
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

      ## single-path search:
      for(j in 1:max.regs){

        ## matrices:
        mXdell <- if(length(delete.adj)==0){NULL}else{mX[,delete.adj]}
        mXndell <- mX[,keep.adj]

        ## estimate model:
        mXadj <- cbind(mXdell,mXndell)
        est <- ols(loge2, mXadj, tol = tol, LAPACK=LAPACK,
          method=2)
        fit <- as.vector(mXadj%*%cbind(est$coefficients))
        ustar <- as.vector(loge2) - fit
        ustar2 <- ustar^2
        mXadj.k <- ncol(mXadj)
        d.f. <- mXadj.n - mXadj.k
        sumUstar2 <- sum(ustar2)
        sigma2 <- sumUstar2/d.f.

        ## make standardised residuals {z_hat}:
        Elnz2 <- -log(mean(exp(ustar)))
        vconstadj <- -Elnz2
        sigma2adj <- exp(fit + vconstadj)
        zhat <- eadj/sqrt(sigma2adj)

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

        ## if diagnostics are ok (i.e. TRUE):
        if(diagnostics.chk){

          #stop if no more deletable regressors:
          if(length(delete.adj)==0){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          #estimate s.e.; compute t-stats. and p-vals.:
          if(vcov.type == "ordinary"){
            coef.var <-as.vector(sigma2*diag(est$xtxinv))
            s.e. <- sqrt(coef.var)
          }
          t.stat <- est$coefficients/s.e.
          p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

          #if any p-value > t.pval:
          if(sum(p.val[1:I(length(delete.adj))] > t.pval) > 0){

            reg.no <- which.max(p.val[1:I(length(delete.adj))])

            #do pet test:
            if(do.pet){
              deleted <- setdiff(delete, delete.adj[-reg.no])
              n.deleted <- length(deleted)
              mR <- NULL #initiate restriction matrix
              for(k in 1:gum.n){
                if(gum.regs[k] %in% deleted){
                  mR <- rbind(mR, c(rep(0,I(k-1)), 1, rep(0, I(gum.n-k) )))
                } #end if(gum.regs[k}..)
              } #end for(k in ..)
              mRestq <- mR%*%cbind(gum.coefs)
              wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=tol)%*%mRestq
              pet.chk <- as.logical(wald.pval < pchisq(wald.stat, n.deleted, lower.tail = FALSE))
            }else{
              pet.chk <- TRUE
            } #end if(do.pet)else..

            #delete regressor if(pet.chk), else move to keep:
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
          } #end if..else.. p-value > t.pval, then delete regressor

        } #end if diagnostics are ok

      } #### end single-path search

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj is already in spec:
      for(l in 1:length(spec)){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break} #stop for(l in..)
      }

      #if spec.adj not in spec (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to spec:
        spec[[length(spec)+1]] <- spec.adj

        #specification results
        n.spec.adj <- length(spec.adj)
        logl <- -eadj.n*log(2*pi)/2 - sum(log(sigma2adj))/2 - sum(eadj^2/sigma2adj)/2
        info.results <- info.criterion(logl, eadj.n,
          n.spec.adj, method=info.method)

        #add terminal to spec:
        spec.results <- rbind(spec.results,
          c(info.results$value, logl, info.results$n,
          info.results$k))
        row.labels <- c(row.labels, paste("spec", length(spec), sep=""))

      } #end if(chk.spec==FALSE)

    } #end multi-path search: for(i in 1:n.paths) loop
  } #end if(n.paths > 0)
} #end if(ar/arch.chk and delete.n>1)

  ## FIND THE BEST MODEL: ##############

  if(!is.null(spec.results)){

    J <- 1:nrow(spec.results)
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

  ## if no search has been undertaken:
  if(is.null(spec.results)){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }

  ## if search has been undertaken:
  if(!is.null(spec.results)){

    ## terminals results:
    if(length(paths)==0){
      out$paths <- NULL
    }else{ out$paths <- paths }
    out$terminals <- spec
    colnames(spec.results) <- col.labels
    where.empty <- which(spec.results[,"k"]==keep.n)
    if(include.empty==FALSE && length(where.empty) > 0){
      row.labels[where.empty] <- paste("spec ", where.empty,
        " (empty):", sep="")
    }
    rownames(spec.results) <- row.labels
    out$terminals.results <- spec.results

    ##if not estimate specific:
    if(!estimate.specific){
      if(best.spec==0 || is.na(best.spec) || length(best.spec)==0 ){
        out$specific.spec <- NULL
      }else{
        specific <- sort(best.spec)
        names(specific) <- object$aux$vXnames[specific]
        out$specific.spec <- specific
      }
    } #end if(!estimate.specific)

    ## if estimate specific:
    if(estimate.specific){

      ## prepare estimation:
      e <- zoo(cbind(eadj), order.by=eadj.index)
      colnames(e) <- "e"
      specific <- sort(best.spec)
      specificadj <- setdiff(specific,1)
      if(length(specificadj)==0){
        vXadj <- NULL
      }else{
        vXadj <- cbind(object$aux$vX[,specificadj])
        colnames(vXadj) <- object$aux$vXnames[specificadj]
        vXadj <- zoo(vXadj, order.by=eadj.index)
      }
      if(is.null(ar.LjungB)){
        ar.LjungB <- object$aux$qstat.options[1]
      }
      if(is.null(arch.LjungB)){
        arch.LjungB <- object$aux$qstat.options[2]
      }

      ## estimate model:
      est <- arx(e, vc=TRUE, vxreg=vXadj,
        zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        tol=object$aux$tol, LAPACK=object$aux$LAPACK,
        verbose=TRUE, plot=FALSE)

      ## rename various stuff:
      est$call <- est$date <- NULL
      est$diagnostics <- est$diagnostics[1:3,]
      where.diagnostics <- which(names(est)=="diagnostics")
      if(length(where.diagnostics)>0){
        names(est)[where.diagnostics] <- "specific.diagnostics"
      }
      est$aux$y.name <- "e"

      ## delete various stuff:
      est$aux$vxreg <- NULL
      est$aux$vxreg.index <- NULL
      est$mean.fit <- NULL

      ## finalise:
      est <- unclass(est)
      names(specific) <- est$aux$vXnames
      out$specific.spec <- specific
      out <- c(out,est)

    } #end if(estimate.specific)

  } #end if(!is.null(spec.results))

  if(length(notes) > 0){ out$notes <- notes }
  out$aux$vXnames.gum <- object$aux$vXnames
  if(is.null(out$aux$vcov.type)){ out$aux$vcov.type <- vcov.type }
  out <- c(list(date=date(), gets.type="getsv"), out)
  out$time.finished <- date()
  class(out) <- "gets"

  if(alarm){ alarm() }
  if(plot){ plot.gets(out) }
  return(out)
}
