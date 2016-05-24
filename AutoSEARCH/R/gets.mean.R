gets.mean <-
function(y, mc=NULL, ar=NULL, ewma=NULL, mx=NULL,
  arch=NULL, asym=NULL, log.ewma=NULL, vx=NULL, keep=NULL, p=2,
  varcov.mat=c("ordinary", "white"), t.pval=0.05, do.pet=TRUE,
  wald.pval=0.05, ar.LjungB=c(2,0.025), arch.LjungB=c(2,0.025),
  tau=2, info.method=c("sc", "aic", "hq"),
  info.resids=c("mean", "standardised"), include.empty=FALSE,
  zero.adj=0.1, vc.adj=TRUE, tol=1e-07, LAPACK=FALSE, max.regs=1000,
  verbose=TRUE, smpl=NULL, alarm=FALSE)
{
### ARGUMENTS ###########

keep.groups=NULL

#match arguments:
info.method <- match.arg(info.method)
info.resids <- match.arg(info.resids)
varcov.mat <- match.arg(varcov.mat)

##zoo:
y <- as.zoo(y)
if(is.null(smpl)){
  iStart <- start(y); iEnd <- end(y)
}else{
  iStart <- smpl[1]; iEnd <- smpl[2]
}
y <- window(as.zoo(y), start=iStart, end=iEnd)
y <- na.trim(y)
zoo.index.y <- index(y)
y <- coredata(y)
y.n <- length(y)
t1 <- zoo.index.y[1]
t2 <- zoo.index.y[y.n]

if(!is.null(mx)){
  mx <- as.zoo(mx)
  mx <- window(mx, start=t1, end=t2)
  mx <- cbind(coredata(mx))
}

if(!is.null(vx)){
  vx <- as.zoo(vx)
  vx <- window(vx, start=t1, end=t2)
  vx <- cbind(coredata(vx))
}

### INITIALISE ##########

out <- list()
out$call <- sys.call()
warnings <- list()
spec <- list()
spec.results <- NULL

#check if variance equation:
if(!is.null(arch) || !is.null(asym) || !is.null(log.ewma) || !is.null(vx)){
  var.spec.chk <- TRUE
}else{var.spec.chk <- FALSE}

## GUM ############################################

##mean regressors:
ewma.mean.chk <- if(is.null(ewma)){0}else{ifelse(is.null(ewma$lag),1,ewma$lag)}
max.ar <- if( is.null(ar)&&is.null(ewma) ){0}else{ max(ar,ewma.mean.chk) }
#max.ar <- if(is.null(ar)){0}else{max(ar)}
mXunadj <- regs.mean.sm(y, mc=mc, ar=ar, ewma=ewma, mx=mx)
y.n <- length(y)
if(is.null(mXunadj)){
  stop("Mean equation empty")
}else{
  mX <- cbind(mXunadj[I(max.ar+1):y.n,])
}
yadj <- y[I(max.ar+1):y.n]

#misc.:
yadj.n <- length(yadj)
keep.n <- length(keep)
gum.n <- ncol(mX)
gum <- 1:gum.n
delete <- setdiff(gum, keep)
delete.n <- length(delete)

#deletable and non-deletable regressors:
if(delete.n > 0){mXdel <- cbind(mX[,delete])}else{mXdel <- NULL}
if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(mX[,keep])}

#estimate GUM:
mXadj <- cbind(mXdel,mXndel)
est <- ols.fit2(yadj, mXadj, tol = tol, LAPACK=LAPACK)
fit <- as.vector(mXadj%*%cbind(est$coefficients))
resids <- as.vector(yadj) - fit
resids.n <- length(resids)
resids2 <- resids^2
mXadj.k <- NCOL(mXadj)
d.f. <- yadj.n - mXadj.k
sumResids2 <- sum(resids2)
sigma2 <- sumResids2/d.f.

#estimate s.e.; compute t-stats. and p-vals.:
if(varcov.mat == "ordinary"){
  varcovmat <- sigma2*est$xtxinv
  coef.var <-as.vector(diag(varcovmat))
  s.e. <- sqrt(coef.var)
}
if(varcov.mat == "white"){
  matResids2 <- matrix(0, yadj.n, yadj.n)
  diag(matResids2) <- resids2
  omega.hat <- t(mXadj)%*%matResids2%*%mXadj
  varcovmat <- est$xtxinv%*%omega.hat%*%est$xtxinv
  coef.var <- as.vector(diag(varcovmat))
  s.e. <- sqrt(coef.var)
}
t.stat <- est$coefficients/s.e.
p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

#make standardised residuals {z_hat}:
if(var.spec.chk){
  resids.n <- length(resids)
  if(!is.null(vx)){
    vx.names <- colnames(vx)
    vx <- cbind(vx[I(yadj.n-resids.n+1):yadj.n,])
    colnames(vx) <- vx.names
  }
  est.var <- sm(resids, mc=NULL, arch=arch, asym=asym,
    log.ewma=log.ewma, vx=vx, p=p, zero.adj=zero.adj,
    vc.adj=vc.adj, tol=tol, LAPACK=LAPACK, verbose=FALSE)
  zhat <- coredata(na.trim(est.var$resids.std))
}else{
  zhat <- resids/sqrt(sigma2)
}

#make diagnostics table:
#if(!is.null(ar.LjungB) || !is.null(arch.LjungB)){
  if(verbose == TRUE){
    diagnostics <- matrix(NA, 2, 3)
    colnames(diagnostics) <- c("Chi^2", "df", "p-val")
    rownames(diagnostics) <- c(paste("Ljung-Box AR(", ar.LjungB[1],
      ")", sep=""), paste("Ljung-Box ARCH(", arch.LjungB[1], ")",
      sep=""))
  }
#}

#Ljung-Box test for serial correlation in {z_t}:
if(!is.null(ar.LjungB)){
  ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
#  ar.LjungBox <- LjungBox(zhat, lag = ar.LjungB[1], zero.mean = TRUE)
  if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.gum.chk <- 0}else{ar.gum.chk <- 1}
  if(verbose == TRUE){
    diagnostics[1,1] <- ar.LjungBox$statistic
    diagnostics[1,2] <- ar.LjungB[1]
    diagnostics[1,3] <- ar.LjungBox$p.value
  }
}else{ar.gum.chk <- 1}

#Ljung-Box test for arch in {z_t^2}
if(!is.null(arch.LjungB)){
  arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
#  arch.LjungBox <- LjungBox(zhat^2, lag = arch.LjungB[1], zero.mean = FALSE)
  if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.gum.chk <- 0}else{arch.gum.chk <- 1}
  if(verbose == TRUE){
    diagnostics[2,1] <- arch.LjungBox$statistic
    diagnostics[2,2] <- arch.LjungB[1]
    diagnostics[2,3] <- arch.LjungBox$p.value
  }
}else{arch.gum.chk <- 1}

#results:
if(verbose == TRUE){
  keep.labels <- c(rep(0,delete.n), rep(1,keep.n))
  if(!is.null(keep.groups)){
    keep.groups.n <- length(keep.groups)
    if(keep.groups.n < 10){divisor <- 10}else{divisor <- 100}
    for(g in 1:length(keep.groups)){
      keep.groups.g <- keep.groups[[g]]
      keep.groups.g.n <- length(keep.groups.g)
      tmp.which <- which(c(delete,keep) %in% keep.groups.g[2:keep.groups.n])
      keep.labels[tmp.which] <- 1 + g/divisor
    }
  }

  results <- as.data.frame(cbind(c(delete,keep), keep.labels,
    est$coefficients, s.e., t.stat, p.val))
  colnames(results) <- c("reg.no", "keep", "coef", "s.e.", "t-stat", "p-val")
  rownames(results) <- colnames(mXunadj[,c(delete,keep)]) #NULL
  out$gum.mean <- results

  if(var.spec.chk){
    out$gum.volatility <- est.var$variance.results
  } #end if(var.spec.chk)

  out$gum.diagnostics <- diagnostics
}

#if GUM passes diagnostic checks:
if((ar.gum.chk*arch.gum.chk) != 0){

  if(info.resids=="standardised" && var.spec.chk==FALSE){
    warnings <- c(warnings, c("Mean residuals used in information criterion"))
  }

  spec[[1]] <- gum

  #specification results
  if(info.resids=="mean" || var.spec.chk==FALSE){
    T <- resids.n
    logl <- I(-1)*T*(log(2*pi) + log(sigma2))/2 - sumResids2/(sigma2*2)
    info.results <- info.criterion(logl, T, mXadj.k, method = info.method)
  }else{
    T <- length(zhat)
    vXadj.k <- NROW(est.var$variance.results["coef"])
    phat <- if(is.null(tau)){gedestp(zhat)}else{tau}
    logl <- gedlogl(zhat, p = phat)
    info.results <- info.criterion(logl, T, I(mXadj.k+vXadj.k),
      method = info.method)
  }
  spec.results <- rbind( c(info.results$value, info.results$n, info.results$k) )
  col.labels <- c(paste("info(", info.method, ")", sep=""), "T", "k")
  row.labels <- c("spec1 (gum)")

  #record data for tests against gum:
  gum.regs <- c(delete, keep)
  gum.coefs <- est$coefficients
  gum.varcovmat <- varcovmat

}else{
  warnings <- c(warnings, c("MGUM does not pass one or more diagnostic checks"))
}

## EMPTY MODEL ####################################

if( ar.gum.chk*arch.gum.chk!=0 && delete.n>0 && include.empty==TRUE ){

  #benchmark model:
  if(keep.n==0){
    resids <- as.vector(yadj)
    resids.n <- length(resids)
    fit <- rep(0, yadj.n)
    mXndel.k <- 0
    resids2 <- resids^2; sumResids2 <- sum(resids2)
    sigma2 <- sumResids2/yadj.n
  }
  if(keep.n>0){
    est <- ols.fit1(yadj, mXndel, tol = tol, LAPACK=LAPACK)
    fit <- as.vector(mXndel%*%cbind(est$coefficients))
    resids <- as.vector(yadj) - fit
    mXndel.k <- NCOL(mXndel)
    resids2 <- resids^2; sumResids2 <- sum(resids2)
    sigma2 <- sumResids2/(yadj.n-mXndel.k)
  }

  #estimate zhat:
  if(var.spec.chk){
    est.var <- sm(resids, mc=NULL, arch=arch, asym=asym,
      log.ewma=log.ewma, vx=vx, zero.adj=zero.adj, vc.adj=vc.adj,
      p=p, tol=tol, LAPACK=LAPACK, verbose=FALSE)
    zhat <- coredata(na.trim(est.var$resids.std))
  }else{
    zhat <- resids/sqrt(sigma2)
  }

  #Ljung-Box test for serial correlation in {z_t}:
  if(!is.null(ar.LjungB)){
    ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
    if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.bench.chk <- 0}else{ar.bench.chk <- 1}
  }else{ar.bench.chk <- 1}

  #Ljung-Box test for arch in {z_t}
  if(!is.null(arch.LjungB)){
    arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.bench.chk <- 0}else{arch.bench.chk <- 1}
  }else{arch.bench.chk <- 1}

  #check if empty model passes diagnostic checks:
  if((ar.bench.chk*arch.bench.chk) != 0){

    spec[[length(spec)+1]] <- if(is.null(keep)){0}else{keep}

    #specification results
    if(info.resids=="mean" || var.spec.chk==FALSE){
      T <- resids.n
      logl <- I(-1)*T*(log(2*pi) + log(sigma2))/2 - sumResids2/(sigma2*2)
      info.results <- info.criterion(logl, T, keep.n, method = info.method)
    }else{
      T <- length(zhat)
      vXadj.k <- NROW(est.var$variance.results["coef"])
      phat <- if(is.null(tau)){gedestp(zhat)}else{tau}
      logl <- gedlogl(zhat, p = phat)
      info.results <- info.criterion(logl, T, I(keep.n+vXadj.k),
        method = info.method)
    }

    #compute wald-test:

    #check pet result:
      spec.results <- rbind(spec.results, c(info.results$value,
        info.results$n, info.results$k))
      row.labels <- c(row.labels, paste("spec", length(spec), " (empty)", sep=""))

  }else{
    warnings <- c(warnings, c("Empty mean model does not pass one or more diagnostic checks"))
  } #end if(empty passes diagnostics==TRUE){..}else{..}
} #end if(include empty model==TRUE)

## MULTI-PATH SEARCH #################

insig.regs <- NULL
paths <- list()
if( ar.gum.chk*arch.gum.chk!=0 && delete.n>1 ){

  #paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  n.paths <- length(insig.regs)

  #if paths = 0:
  if(n.paths == 0){
    warnings <- c(warnings, c("No insignificant regressors in MGUM"))
  }

  #if paths > 0:
  if(n.paths > 0){

    #paths:
    for(i in 1:n.paths){

      #print path if verbose:
      if(verbose){
        print(paste("Searching path no. ", i, " out of ", n.paths, sep=""),
        quote=FALSE)
      }

      #prepare single-path search:
      path <- insig.regs[i]
      delete.adj <- setdiff(delete, insig.regs[i])
      keep.adj <- as.numeric(keep)

      #single-path search:
      for(j in 1:max.regs){

        #matrices:
        mXdell <- if(length(delete.adj)==0){NULL}else{mX[,delete.adj]}
        mXndell <- if(is.null(keep.adj)){NULL}else{mX[,keep.adj]}

        #estimate model:
        mXadj <- cbind(mXdell,mXndell)
        if(!is.null(mXadj)){
          est <- ols.fit2(yadj, mXadj, tol = tol, LAPACK=LAPACK)
          fit <- as.vector(mXadj%*%cbind(est$coefficients))
        }else{
          fit <- rep(0, yadj.n)
        }
        resids <- as.vector(yadj) - fit
        resids2 <- resids^2
        mXadj.k <- ncol(mXadj)
        d.f. <- yadj.n - mXadj.k
        sumResids2 <- sum(resids2)
        sigma2 <- sumResids2/d.f.

        #make standardised residuals {z_hat}:
        if(var.spec.chk){
          est.var <- sm(resids, mc=NULL, arch=arch, asym=asym,
            log.ewma=log.ewma, vx=vx, zero.adj=zero.adj, vc.adj=vc.adj,
            p=p, tol=tol, LAPACK=LAPACK, verbose=FALSE)
          zhat <- coredata(na.trim(est.var$resids.std))
        }else{
          zhat <- resids/sqrt(sigma2)
        }

        #Ljung-Box test for serial correlation in {z_t}:
        if(!is.null(ar.LjungB)){
          ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
          if(ar.LjungBox$p.value <= ar.LjungB[2]){ar.chk <- 0}else{ar.chk <- 1}
        }else{ar.chk <- 1}

        #Ljung-Box test for arch in {z_t}
        if(!is.null(arch.LjungB)){
          arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
          if(arch.LjungBox$p.value <= arch.LjungB[2]){arch.chk <- 0}else{arch.chk <- 1}
        }else{arch.chk <- 1}

        #if either ar.chk or arch.chk is equal to 0, then move
        #path[length(path)] to keep.adj
        if((ar.chk*arch.chk) == 0){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*I(-1))
          next #next j
        }

        #check keep.groups:
        if(!is.null(keep.groups)){
          skip.chk <- FALSE
          for(g in 1:length(keep.groups)){
            group.g.n <- length(keep.groups[[g]])
            group.g <- keep.groups[[g]][2:group.g.n]
            vars.in.delete <- intersect(group.g, delete.adj)
            vars.in.delete.n <- length(vars.in.delete)
            vars.in.keep <- intersect(group.g, keep.adj)
            vars.in.keep.n <- length(vars.in.keep)
            if(vars.in.delete.n + vars.in.keep.n == keep.groups[[g]][1]){
              if(vars.in.delete.n > 0){
                keep.adj <- union(vars.in.delete, keep.adj)
                delete.adj <- setdiff(delete.adj, vars.in.delete)
                path <- union(path, vars.in.delete*I(-1))
                skip.chk <- TRUE
                break  ##break for loop with index g
              }
            }
          }
          if(skip.chk == TRUE){next}  ##next j
        } #end check keep.groups

        #if ar.chk*arch.chk==1:
        if(ar.chk*arch.chk == 1){

          #stop if no more deletable regressors:
          if(length(delete.adj)==0){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          if(!is.null(mXadj)){

            #estimate s.e.; compute t-stats. and p-vals.:
            if(varcov.mat == "ordinary"){
              coef.var <-as.vector(sigma2*diag(est$xtxinv))
              s.e. <- sqrt(coef.var)
            } #end "ordinary"

            if(varcov.mat == "white"){
              matResids2 <- matrix(0, yadj.n, yadj.n)
              diag(matResids2) <- resids2
              omega.hat <- t(mXadj)%*%matResids2%*%mXadj
              varcovmat <- est$xtxinv%*%omega.hat%*%est$xtxinv
              coef.var <- as.vector(diag(varcovmat))
              s.e. <- sqrt(coef.var)
            } #end "white"

            #t-tests:
            t.stat <- est$coefficients/s.e.
            p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

          } #end if(!is.null(mXadj))

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
            } #end if(do.pet)
            #if pet=TRUE then delete regressor, else move to keep:
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

        } #end if ar.chk*arch.chk==1

      } #### end single-path search

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj is already in spec:
      n.spec <- length(spec)
      if(length(spec.adj)==0){spec.adj <- 0}
      for(l in 1:n.spec){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break}
      }

      #if spec.adj not in spec (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to spec:
        spec[[length(spec)+1]] <- spec.adj

        #specification results
        n.spec.adj <- length(spec.adj)
        if(info.resids=="mean" || var.spec.chk==FALSE){
          T <- resids.n
          logl <- I(-1)*T*(log(2*pi) + log(sigma2))/2 - sumResids2/(sigma2*2)
          info.results <- info.criterion(logl, T, n.spec.adj, method = info.method)
        }else{
          T <- length(zhat)
          phat <- if(is.null(tau)){gedestp(zhat)}else{tau}
          logl <- gedlogl(zhat, p = phat)
          info.results <- info.criterion(logl, T, I(n.spec.adj+vXadj.k),
            method=info.method)
        } #end specification results

        #add terminal to spec:
        spec.results <- rbind(spec.results, c(info.results$value,
          info.results$n, info.results$k))
        row.labels <- c(row.labels, paste("spec", length(spec), sep=""))

      } #end if(chk.spec==FALSE)
    } #end multi-path search: for(i in 1:n.paths) loop
  } #end if paths > 0
} #end if(ar/arch.chk and delete.n>1)

## FIND THE BEST MODEL ########################

if(verbose){
  if(!is.null(spec.results)){

    #if keep.groups != NULL:
    if(!is.null(keep.groups)){
      if((ar.bench.chk*arch.bench.chk) == 1){
        J <- I(1:NROW(spec.results))[-2]
        models <- cbind(J, spec.results[-2,])
        colnames(models) <- NULL
        warnings <- c(warnings, "keep.groups not NULL: Empty model (spec2) not included in selection")
      }else{
        J <- 1:nrow(spec.results)
        models <- cbind(J, spec.results)
        colnames(models) <- NULL
      }
    } #end if(!is.null(keep.groups))

    #if keep.groups == NULL
    if(is.null(keep.groups)){
      J <- 1:nrow(spec.results)
      models <- cbind(J, spec.results)
      colnames(models) <- NULL
    } #end if(is.null(keep.groups))

    #find best model:
    where <- which.min(models[,2])
    best.spec <- spec[[where]] #winner

    #check for several minimums:
    min.models <- min(models[,2])
    wheres <- which(models[,2]==min.models)
    if(length(wheres)>1){warnings <- c(warnings, "Several terminal specifications attain the minimum information criterion")}

    #best specification:
    #best.spec <- spec[[winner]]

  } #end if(!is.null(spec.results))
} #end if(verbose)

## ESTIMATE SPECIFIC MODEL #############

specific.mean <- NULL
specific.diagnostics <- NULL
if(verbose){
  if(!is.null(spec.results)){

    if(best.spec[1]==0){
      resids <- yadj
      specific.mean <- "spec2 (empty)"
      #R-squared:
      if(verbose){ Rsquared <- 0 }
    }else{

      #estimate specific:
      specific <- sort(best.spec)  ##specific
      mXadj <- cbind(mX[,specific])
      est <- ols.fit2(yadj, mXadj, tol = tol, LAPACK=LAPACK)
      fit <- as.vector(mXadj%*%cbind(est$coefficients))
      resids <- as.vector(yadj) - fit
      resids2 <- resids^2
      mXadj.k <- ncol(mXadj)
      d.f. <- yadj.n - mXadj.k
      sumResids2 <- sum(resids2)
      sigma2 <- sumResids2/d.f.
      #R-squared:
      if(verbose){
        TSS <- sum( (yadj - mean(yadj))^2 )
        RSS <- sum( (resids - mean(resids))^2 )
        Rsquared <- 1 - RSS/TSS
      }

      #estimate s.e.; compute t-stats. and p-vals.:
      if(varcov.mat == "ordinary"){
        coef.var <-as.vector(sigma2*diag(est$xtxinv))
        s.e. <- sqrt(coef.var)
      }
      if(varcov.mat == "white"){
        matResids2 <- matrix(0, yadj.n, yadj.n)
        diag(matResids2) <- resids2
        omega.hat <- t(mXadj)%*%matResids2%*%mXadj
        varcovmat <- est$xtxinv%*%omega.hat%*%est$xtxinv
        coef.var <- as.vector(diag(varcovmat))
        s.e. <- sqrt(coef.var)
      }
      t.stat <- est$coefficients/s.e.
      p.val <- pt(abs(t.stat), d.f., lower.tail=FALSE)*2

      #make results table:
      specific.mean <- as.data.frame(cbind(specific, est$coefficients, s.e., t.stat, p.val))
      colnames(specific.mean) <- c("reg", "coef", "s.e.", "t-stat", "p-val")
      rownames(specific.mean) <- colnames(mXunadj)[specific.mean[,1]] #NULL

    } #end if(best.spec[1]==0)else{..}

    #make standardised residuals {z_hat}:
    if(var.spec.chk){
      est.var <- sm(resids, mc=NULL, arch=arch, asym=asym,
        log.ewma=log.ewma, vx=vx, p=p, zero.adj=zero.adj,
        vc.adj=vc.adj, tol=tol, LAPACK=LAPACK, verbose=FALSE)
      specific.variance <- cbind(est.var$variance.results)
      zhat <- est.var$resids.std
      zhat <- zoo(c(rep(NA, I(y.n-length(zhat))), zhat), order.by=zoo.index.y)
      resids <- zoo(c(rep(NA, I(y.n-length(resids))), resids), order.by=zoo.index.y)
      #out <- c(list(resids=resids, resids.std=zhat), out)
      if(p==2){
        out <- c(list(resids=resids, resids.std=zhat, Elogzp=est.var$Elogzp), out)
      }else{
        out <- c(list(resids=resids, resids.std=zhat, Elogzp=est.var$Elogzp,
          Elogzstarp=est.var$Elogzstarp, logEzp=est.var$logEzp), out)
      }

    }else{
      zhat <- resids/sqrt(sigma2)
      zhat <- zoo(c(rep(NA, I(y.n-length(zhat))), zhat), order.by=zoo.index.y)
      resids <- zoo(c(rep(NA, I(y.n-length(resids))), resids), order.by=zoo.index.y)
      out <- c(list(resids=resids, resids.std=zhat), out)
    } #end if(var.spec.chk){..}else{..}

    #make diagnostics table:
    if(!is.null(ar.LjungB) || !is.null(arch.LjungB)){
      #specific.diagnostics <- matrix(NA, 2, 3)
      specific.diagnostics <- matrix(NA, 3, 3)
      colnames(specific.diagnostics) <- c("Chi^2", "df", "p-val")
      rownames(specific.diagnostics) <- c(paste("Ljung-Box AR(", ar.LjungB[1],
        ")", sep=""), paste("Ljung-Box ARCH(", arch.LjungB[1], ")",
        sep=""), "R-squared")
      if(!is.null(ar.LjungB)){
        ar.LjungBox <- Box.test(zhat, lag = ar.LjungB[1], type="L")
        specific.diagnostics[1,1] <- ar.LjungBox$statistic
        specific.diagnostics[1,2] <- ar.LjungB[1]
        specific.diagnostics[1,3] <- ar.LjungBox$p.value
      }
      if(!is.null(arch.LjungB)){
        arch.LjungBox <- Box.test(zhat^2, lag = arch.LjungB[1], type="L")
        specific.diagnostics[2,1] <- arch.LjungBox$statistic
        specific.diagnostics[2,2] <- arch.LjungB[1]
        specific.diagnostics[2,3] <- arch.LjungBox$p.value
      }
      ##R-squared:
      specific.diagnostics[3,1] <- Rsquared

    } #end make diagnostics table
  } #end if(!is.null(spec.results))
} #end if(verbose)

## CREATE OUT LIST ################################

out$keep <- keep
out$insigs.in.gum <- insig.regs
#out$smpl <- c(t1,t2)

if((ar.gum.chk*arch.gum.chk) != 0){
  if(verbose){ if(length(paths)==0){out$paths <- NULL}else{out$paths <- paths} }
  out$terminals <- spec
  colnames(spec.results) <- col.labels
  rownames(spec.results) <- row.labels
  out$terminals.results <- spec.results
  if(verbose){
    out$specific.mean <- specific.mean
    if(var.spec.chk){out$specific.volatility <- specific.variance}
    out$specific.diagnostics <- specific.diagnostics
  }
} #end if((ar.gum.chk*arch.gum.chk) != 0)

if(length(warnings) > 0){
  out$warnings <- warnings
}
if(alarm){alarm()}

return(out)
}
