soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}
 
mean.na <- function(vec){
  return(mean(vec[!is.na(vec)]))
}

l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}

safesvd <- function(x){
  i <- 1
  out <- try(svd(x), silent=TRUE)
  while(i<10 && class(out)=="try-error"){
    out <- try(svd(x), silent=TRUE)
    i <- i+1
  }
  if(class(out)=="try-error") out <- svd(matrix(rnorm(nrow(x)*ncol(x)), ncol=ncol(x)))
  return(out)
}

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}



SMD <- function(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg){
  # This gets a single factor. Do MultiSMD to get multiple factors.
  nas <- is.na(x)
  v.init <- v
  xoo <- x
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  oldv <- rnorm(ncol(x))
  for(iter in 1:niter){
    if(sum(abs(oldv-v))>1e-7){
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      lamu <- BinarySearch(argu,sumabsu)
      su <- soft(argu,lamu)
      u <- matrix(su/l2n(su),ncol=1)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      lamv <- BinarySearch(argv, sumabsv)
      sv <- soft(argv,lamv)
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    } 
  }
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  return(list(d=d, u=u, v=v, v.init=v.init))
}

SMDOrth <- function(x, us, sumabsv=NULL, niter=20, trace=TRUE,v, vpos, vneg){
  # Gets the next u for sparse PCA, using Trevor's method of requiring this new u to be orthog
  # to all previous us (contained in us)
  nas <- is.na(x)
  v.init <- v
  xoo <- x
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  u <- rnorm(nrow(x))
  oldu <- rnorm(nrow(x))
  oldv <- rnorm(ncol(x))
  for(iter in 1:niter){
    if(sum(abs(oldu-u))>1e-6 || sum(abs(oldv-v))>1e-6){
      oldu <- u
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      numer <- lsfit(y=argu, x=us, intercept=FALSE)$res
      u <- numer/l2n(numer)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      lamv <- BinarySearch(argv, sumabsv)
      sv <- soft(argv,lamv)
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    } 
  }
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  return(list(d=d,u=u,v=v, v.init=v.init))
}


CheckPMDV <-  function(v,x,K){
  if(!is.null(v) && is.matrix(v) && ncol(v)>=K){
    v <- matrix(v[,1:K], ncol=K)
  } else if(ncol(x)>nrow(x)){
    x[is.na(x)] <- mean.na(x)
    v <- matrix(t(x)%*%(safesvd(x%*%t(x))$v[,1:K]),ncol=K)
    if(sum(is.na(v))>0) v <- matrix(safesvd(x)$v[,1:K], ncol=K)
    v <- sweep(v,2,apply(v, 2, l2n), "/")
    if(sum(is.na(v))>0) stop("some are NA")
  } else if (ncol(x)<=nrow(x)){
    x[is.na(x)] <- mean.na(x)
    v <- matrix(safesvd(t(x)%*%x)$v[,1:K],ncol=K)
  }
  return(v)
}

SPC <- function(x, sumabsv=4, niter=20, K=1, orth=FALSE, trace=TRUE, v=NULL, center=TRUE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE){
  if(vpos&&vneg) stop("Cannot constrain elements to be positive AND negative.")
  out <- PMDL1L1(x,sumabsu=sqrt(nrow(x)), sumabsv=sumabsv, niter=niter,K=K,orth=orth,trace=trace,v=v,center=center,cnames=cnames, upos=FALSE, uneg=FALSE, vpos=vpos, vneg=vneg)
  if(compute.pve){
     # Calculate percent variance explained, using definition on page 7 of Shen and Huang (2008) Journal of Multivariate Analysis vol 99: 1015-1034
    v <- matrix(out$v, ncol=K)
    ve <- NULL # variance explained
    xfill <- x
    if(center) xfill <- x-mean.na(x)
    xfill[is.na(x)] <- mean.na(xfill)
    for(k in 1:K){
      vk <- matrix(v[,1:k], ncol=k)
      xk <- xfill%*%vk%*%solve(t(vk)%*%vk)%*%t(vk)
      svdxk <- svd(xk)
      ve <- c(ve, sum(svdxk$d^2))    
    }
    pve <- ve/sum(svd(xfill)$d^2) # proportion of variance explained
    out$prop.var.explained <- pve
  }
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "SPC"
  return(out)
}

print.SPC <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero v's: ", apply(x$v!=0, 2, sum), "\n")
  cat("Sumabsv used: ", x$sumabsv, "\n")
  cat("Proportion of variance explained by first k components: \n")
  tab <- cbind(paste(paste("First ", sep="", 1:length(x$prop.var.explained)), sep="", " Components"), round(x$prop.var.explained, 3))
  dimnames(tab) <- list(1:length(x$prop.var.explained), c("", "Prop. Var. Explained"))
  print(tab, quote=FALSE, sep="\t")
  if(!is.null(x$meanx)) cat("Subtracted out mean of x: ", x$meanx, "\n")
  if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      v <- x$v[,k]
      if(is.null(x$cnames)) x$cnames <- 1:length(v)
      cat(fill=T)
      vs <- cbind(x$cnames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Feature Name", "Feature Weight"))
      print(vs, quote=FALSE, sep="\t")
    }
  }
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}

PMD <- function(x, type=c("standard", "ordered"), sumabs=.4, sumabsu=5, sumabsv=NULL, lambda=NULL, niter=20, K=1, v=NULL, trace=TRUE, center=TRUE, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE){
  if(upos&&uneg || vpos&&vneg) stop("Cannot contrain elements to be both positive and negative!")
  if(type=="ordered" && (vpos||vneg)) stop("Cannot constrain signs of v if type is ordered.")
  call <- match.call()
  type <- match.arg(type)
  if(type=="standard"){
    out <- PMDL1L1(x, sumabs=sumabs, sumabsu=sumabsu, sumabsv=sumabsv, niter=niter, K=K, v=v, trace=trace, center=center, rnames=rnames, cnames=cnames, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    class(out) <- "PMDL1L1"
  } else if(type=="ordered"){
    out <- PMDL1FL(x,K=K,sumabsu=sumabsu,lambda=lambda,chrom=chrom,niter=niter, v=v, trace=trace, center=center, rnames=rnames, cnames=cnames, upos=upos, uneg=uneg)
    class(out) <- "PMDL1FL"
  }
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  out$call <- call
  return(out)
}

PMD.cv <- function(x, type=c("standard", "ordered"), sumabss=seq(0.1,0.7,len=10), sumabsus=NULL, lambda=NULL, nfolds=5, niter=5, v=NULL, chrom=NULL, nuc=NULL, trace=TRUE, center=TRUE, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE){
  if(upos&&uneg || vpos&&vneg) stop("Cannot contrain elements to be both positive and negative!")
  if(type=="ordered" && (vpos||vneg)) stop("Cannot constrain signs of v if type is ordered.")
  call <- match.call()
  type <- match.arg(type) 
  if(type=="standard"){
    out <- PMDL1L1.cv(x, sumabss=sumabss, nfolds=nfolds, niter=niter, v=v, trace=trace, center=center, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    class(out) <- "PMDL1L1.cv"
  } else if(type=="ordered"){
    out <- PMDL1FL.cv(x, nfolds=nfolds, niter=niter, lambda=lambda, sumabsus=sumabsus, chrom=chrom, v=v, trace=trace, nuc=nuc, center=center, upos=upos, uneg=uneg)
    class(out) <- "PMDL1FL.cv"
  }
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  out$call <- call
  return(out)
}



PMDL1L1 <- function(x,sumabs=.4,sumabsu=NULL,sumabsv=NULL,niter=20,K=1,v=NULL, trace=TRUE, orth=FALSE, center=TRUE, rnames=NULL, cnames=NULL, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg){
  if(center){
    meanx <- mean.na(x)
    x <- x-meanx
  } else {
    meanx <- NULL
  }
  if(orth){
    if(is.null(sumabsu) || sumabsu < sqrt(nrow(x))){
      orth <- FALSE
      warning("Orth option ignored because sparse PCA results only when sumabsu equals sqrt(nrow(x))")
    }
  }
#  if((is.null(sumabsu) && !is.null(sumabsv)) || (is.null(sumabsv) && !is.null(sumabsu))) warning("Sumabsu and sumabsv BOTH must be input when type=standard. Since only one was given, it was ignored and sumabs was used instead.")
  if(is.null(sumabsu) || is.null(sumabsv)){
    sumabsu <- sqrt(nrow(x))*sumabs
    sumabsv <- sqrt(ncol(x))*sumabs
  }
  call <-  match.call()
  if(trace && abs(mean.na(x)) > 1e-15) warning("PMDL1L1 was run without first subtracting out the mean of x.")
  if(!is.null(sumabsu) && (sumabsu<1 || sumabsu>sqrt(nrow(x)))) stop("sumabsu must be between 1 and sqrt(n)")
  if(!is.null(sumabsv) && (sumabsv<1 || sumabsv>sqrt(ncol(x)))) stop("sumabsv must be between 1 and sqrt(p)")
  v <- CheckPMDV(v,x,K)
  if(K>1 && !orth) out <- (MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  if(K>1 && orth) out <- MultiSMDOrth(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v,  vpos=vpos, vneg=vneg)
  if(K==1) out <- SMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
  obj <- (list(u=out$u,v=out$v, d=out$d, v.init=out$v.init, call=call, meanx=meanx,sumabsu=sumabsu, sumabsv=sumabsv, rnames=rnames, cnames=cnames, K=K, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  class(obj) <- "PMDL1L1"
  return(obj)
}

print.PMDL1L1 <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero u's: ", apply(x$u!=0, 2, sum), "\n")
  cat("Num non-zero v's: ", apply(x$v!=0, 2, sum), "\n")
  cat("Sumabsu used: ", x$sumabsu, "\n")
  cat("Sumabsv used: ", x$sumabsv, "\n")
  if(!is.null(x$meanx)) cat("Subtracted out mean of x: ", x$meanx, "\n")
  if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      u <- x$u[,k]
      v <- x$v[,k]
      if(is.null(x$rnames)) x$rnames <- 1:length(u)
      if(is.null(x$cnames)) x$cnames <- 1:length(v)
      cat(fill=T)
      us <- cbind(x$rnames[u!=0], round(u[u!=0],3))
      dimnames(us) <- list(1:sum(u!=0), c("Row Feature Name", "Row Feature Weight"))
      vs <- cbind(x$cnames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Column Feature Name", "Column Feature Weight"))
      print(us, quote=FALSE, sep="\t")
      cat(fill=T)
      print(vs, quote=FALSE, sep="\t")
    }
  }
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}

print.PMDL1FL <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero u's: ", apply(x$u!=0, 2, sum), "\n")
  cat("Num non-zero v's: ", apply(x$v!=0, 2, sum), "\n")
  cat("Sumabsu used: ", x$sumabsu, "\n")
  cat("Lambda used: ", x$lambda, "\n")
  if(!is.null(x$meanx)) cat("Subtracted out mean of x: ", x$meanx, "\n")
    if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      u <- x$u[,k]
      v <- x$v[,k]
      if(is.null(x$rnames)) x$rnames <- 1:length(u)
      if(is.null(x$cnames)) x$cnames <- 1:length(v)
      cat(fill=T)
      us <- cbind(x$rnames[u!=0], round(u[u!=0],3))
      dimnames(us) <- list(1:sum(u!=0), c("Row Feature Name", "Row Feature Weight"))
      vs <- cbind(x$cnames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Column Feature Name", "Column Feature Weight"))
      print(us, quote=FALSE, sep="\t")
      cat(fill=T)
      print(vs, quote=FALSE, sep="\t")
    }
  }
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
}
 
MultiSMDOrth <- function(x,sumabsu,sumabsv,niter,K, trace, v, vpos, vneg){
  if(sumabsu < sqrt(nrow(x))){
    warning("sumabsu was less than sqrt(nrow(x)), so the orthogonal option was not implemented.")
    return(MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, vpos=vpos, vneg=vneg))
  }
  nas <- is.na(x)
  v.init <- v
  xuse <- x
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  for(k in 1:K){
    if(k==1) out <- SMD(xuse,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=FALSE, uneg=FALSE,  vpos=vpos, vneg=vneg)
    if(k>1) out <- SMDOrth(xuse,us[,1:(k-1)],sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace,  vpos=vpos, vneg=vneg)
    us[,k] <- out$u 
    vs[,k] <- out$v 
    ds[k] <- out$d
  }
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}

MultiSMD <- function(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg){
  nas <- is.na(x)
  v.init <- v
  xuse <- x
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  for(k in 1:K){
    out <- SMD(xuse, sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    us[,k] <- out$u 
    vs[,k] <- out$v 
    ds[k] <- out$d
    res <- xuse - out$d*out$u%*%t(out$v) 
    xuse[!nas] <- res[!nas] # [!nas] is new on July 24 2009
  }
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}

PMDL1L1.cv <- function(x, sumabss=seq(0.1,0.7,len=10), nfolds=5, niter=5, v=NULL, trace=TRUE, center=TRUE, upos, uneg, vpos, vneg){
  if(nfolds<2) stop("Must run at least 2 cross-validation folds.")
  percentRemove <- min(.25, 1/nfolds)
  call <- match.call()
  if(max(sumabss)>1 || min(sumabss)<0) stop("sumabs must be between 0 and 1")
  # PMD only cross-validates the single-factor model.
  # and it requires just one sparsity  measure.
  xfill <- x
  missing <- is.na(x)
  v <- CheckPMDV(v,x,K=1)
  errs <- matrix(NA, nrow=nfolds, ncol=length(sumabss))
  nonzerous <- matrix(NA, nrow=nfolds, ncol=length(sumabss))
  nonzerovs <- matrix(NA, nrow=nfolds, ncol=length(sumabss))
  rands <- matrix(runif(nrow(x)*ncol(x)),ncol=ncol(x))
  for(i in 1:nfolds){
    if(trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i-1)*percentRemove < rands)&(rands < i*percentRemove)
    xrm <- x
    xrm[rm] <- NA 
    for(j in 1:length(sumabss)){
      out <- PMDL1L1(xrm, sumabs=sumabss[j], niter=niter, v=v, trace=FALSE, center=center,K=1, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
      xhat <- as.numeric(out$d)*out$u%*%t(out$v)
      errs[i,j] <- sum(((xhat-x)[rm & !missing])^2)
      nonzerous[i,j] <- sum(out$u!=0)
      nonzerovs[i,j] <- sum(out$v!=0)
    }
  }
  if(trace) cat(fill=TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerous.mean <- apply(nonzerous,2,mean)
  nonzerovs.mean <- apply(nonzerovs,2,mean)
  bestsumabs <- sumabss[which.min(err.means)]
  object <- (list(cv=err.means, cv.error=err.sds,bestsumabs=bestsumabs, nonzerous=nonzerous.mean, nonzerovs=nonzerovs.mean, v.init=v, call=call, sumabss=sumabss, nfolds=nfolds))
  class(object) <- "PMDL1L1.cv"
  return(object) 
}

SPC.cv <- function(x, sumabsvs=seq(1.2,5,len=10), nfolds=5, niter=5, v=NULL, trace=TRUE, orth=FALSE,  center=TRUE, vpos=FALSE, vneg=FALSE){
  if(vpos&&vneg) stop("Cannot constrain elements of v to be both positive and negative.")
  if(nfolds<2) stop("Must run at least 2 cross-validation folds.")
  percentRemove <- min(.25, 1/nfolds)
  call <- match.call()
  if(max(sumabsvs)>sqrt(ncol(x)) || min(sumabsvs)<1) stop("sumabs must be between 1 and sqrt(ncol(x))") 
  xfill <- x
  missing <- is.na(x)
  v <- CheckPMDV(v,x,K=1)
  errs <- matrix(NA, nrow=nfolds, ncol=length(sumabsvs))
  nonzerovs <- matrix(NA, nrow=nfolds, ncol=length(sumabsvs))
  rands <- matrix(runif(nrow(x)*ncol(x)),ncol=ncol(x))
  for(i in 1:nfolds){
    if(trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i-1)*percentRemove < rands)&(rands < i*percentRemove)
    xrm <- x
    xrm[rm] <- NA 
    for(j in 1:length(sumabsvs)){
      out <- SPC(xrm, sumabsv=sumabsvs[j], orth=orth, niter=niter, v=v, trace=FALSE, center=center,K=1, vpos=vpos, vneg=vneg)
      xhat <- as.numeric(out$d)*out$u%*%t(out$v)
      errs[i,j] <- sum(((xhat-x)[rm & !missing])^2)
      nonzerovs[i,j] <- sum(out$v!=0)
    }
  }
  if(trace) cat(fill=TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerovs.mean <- apply(nonzerovs,2,mean)
  bestsumabsv <- sumabsvs[which.min(err.means)]
  bestsumabsv1se <- sumabsvs[min(which(err.means < min(err.means) + err.sds[which.min(err.means)]))]
  object <- (list(cv=err.means, cv.error=err.sds,bestsumabsv=bestsumabsv,nonzerovs=nonzerovs.mean, v.init=v, call=call, sumabsvs=sumabsvs, nfolds=nfolds, bestsumabsv1se=bestsumabsv1se, vpos=vpos, vneg=vneg))
  class(object) <- "SPC.cv"
  return(object) 
}

plot.SPC.cv <- function(x,...){
  sumabsvs <- x$sumabsvs
  err.means <- x$cv
  err.sds <- x$cv.error
  plot(sumabsvs, err.means, ylim=range(c(err.means+err.sds, err.means-err.sds)), main="CV Error for SPC", xlab="Value of sumabsv", ylab="Average SSE", type="l")
  points(sumabsvs, err.means)
  lines(sumabsvs, err.means+err.sds, lty="dashed")
  lines(sumabsvs, err.means-err.sds, lty="dashed")
}



print.SPC.cv <- function(x,...){
  cat("Call:\n")
  dput(x$call)
  cat("\n \n")
  cat('Cross-validation errors: \n')
  mat <- round(cbind(round(x$sumabsvs, 3), x$cv, x$cv.error, x$nonzerovs),3)
  dimnames(mat) <- list(1:length(x$sumabsvs), c("Sumabsvs", "CV Error", "CV S.E.", "# non-zero v's"))
  print(mat, quote=F)
  cat('\n Best sumabsv value (lowest CV error): ', x$bestsumabsv, "\n")
  cat('\n Smallest sumabsv value that has CV error within 1 SE of best CV error: ', x$bestsumabsv1se, "\n")
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}

print.PMDL1L1.cv <- function(x,...){
  cat("Call:\n")
  dput(x$call)
  cat("\n \n")
  cat('Cross-validation errors: \n')
  mat <- round(cbind(round(x$sumabss, 3), x$cv, x$cv.error, x$nonzerous, x$nonzerovs),3)
  dimnames(mat) <- list(1:length(x$sumabss), c("Sumabss", "CV Error", "CV S.E.", "# non-zero u's", "# non-zero v's"))
  print(mat, quote=F)
  cat('\n Best sumabs value (lowest CV error): ', x$bestsumabs, "\n")
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
  if(x$vpos) cat("Elements of v constrained to be positive.", fill=TRUE)
  if(x$vneg) cat("Elements of v constrained to be negative.", fill=TRUE)
}

plot.PMDL1L1.cv <- function(x,...){
  sumabss <- x$sumabss
  err.means <- x$cv
  err.sds <- x$cv.error
  plot(sumabss, err.means, ylim=range(c(err.means+err.sds, err.means-err.sds)), main="CV Error for PMD(L1,L1)", xlab="Value of sumabs", ylab="Average SSE", type="l")
  points(sumabss, err.means)
  lines(sumabss, err.means+err.sds, lty="dashed")
  lines(sumabss, err.means-err.sds, lty="dashed")
}


PMDL1FL.cv <- function(x, nfolds=5, niter=5, lambda=NULL, sumabsus=NULL, chrom=NULL, v=NULL, trace=TRUE, nuc=NULL, center=TRUE, upos, uneg){
  if(is.null(sumabsus)) sumabsus <- seq(1, .7*sqrt(nrow(x)),len=10)
  percentRemove <- min(.25, 1/nfolds)
  call <- match.call()
  if(is.null(chrom)) chrom <- rep(1, ncol(x))
  v <- CheckPMDV(v,x,K=1)
  if(is.null(lambda)) lambda <- ChooseLambda1Lambda2(as.numeric(v[,1]))
  errs <- matrix(NA, nrow=nfolds, ncol=length(sumabsus))
  nonzerous <- matrix(NA, nrow=nfolds, ncol=length(sumabsus))
  nonzerovs <- matrix(NA, nrow=nfolds, ncol=length(sumabsus))
  rands <- matrix(runif(nrow(x)*ncol(x)), nrow=nrow(x))
  missing <- is.na(x)
  for(i in 1:nfolds){
    if(trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i-1)*percentRemove < rands)&(rands < i*percentRemove)
    xrm <- x
    xrm[rm] <- NA
    for(j in 1:length(sumabsus)){
      if(trace) cat(j,fill=F)
      out <- PMDL1FL(xrm, K=1, lambda=lambda,sumabsu=sumabsus[j], niter=niter, chrom=chrom, v=v, trace=FALSE, center=center, upos=upos, uneg=uneg)
      xhat <- as.numeric(out$d)*out$u%*%t(out$v)
      errs[i,j] <- sum(((xhat-x)[rm & !missing])^2)
      nonzerous[i,j] <- sum(out$u!=0)
      nonzerovs[i,j] <- sum(out$v!=0)
    }
  }
  if(trace) cat(fill=TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerous.mean <- apply(nonzerous,2,mean)
  nonzerovs.mean <- apply(nonzerovs,2,mean)
  bestind <- which.min(err.means)
  bestsumabsu <- sumabsus[bestind]
  object <- (list(cv=err.means, lambda=lambda,cv.error=err.sds,bestsumabsu=bestsumabsu, nonzerous=nonzerous.mean, nonzerovs=nonzerovs.mean, v.init=v, call=call, sumabsus=sumabsus, nfolds=nfolds,x=x,chrom=chrom,nuc=nuc, niter=niter))
  class(object) <- "PMDL1FL.cv"
  return(object) 
}

plot.PMDL1FL.cv <- function(x,...){
  mat <- x$x
  sumabsus <- x$sumabsus
  err.means <- x$cv
  err.sds <- x$cv.error
  chrom <- x$chrom
  niter <- x$niter
  nuc <- x$nuc
  lambda <- x$lambda
  v <- x$v.init
  bestsumabsu <- x$bestsumabsu
  par(mfrow=c(2,1))
  plot(sumabsus, err.means, ylim=range(c(err.means+err.sds, err.means-err.sds)), main="CV Error for PMD(L1, FL)", xlab="Value of sumabsu", ylab="Average SSE", type="l")
  points(sumabsus, err.means)
  lines(sumabsus, err.means+err.sds, lty="dashed")
  lines(sumabsus, err.means-err.sds, lty="dashed")
  out <- PMDL1FL(mat, sumabsu=bestsumabsu,  K=1, chrom=chrom, niter=niter, v=v, trace=FALSE, lambda=lambda, upos=x$upos, uneg=x$uneg)
  PlotCGH(out$v, main="Best V", chrom=chrom, nuc=nuc)
}

print.PMDL1FL.cv <- function(x,...){
  cat("Call:\n")
  dput(x$call)
  cat("\n \n")
  cat('Cross-validation errors: \n')
  mat <- round(cbind(round(x$sumabsus, 3), x$cv, x$cv.error, x$nonzerous, x$nonzerovs),3)
  dimnames(mat) <- list(1:length(x$sumabsus), c("Sumabsus", "CV Error", "CV S.E.", "# non-zero u's", "# non-zero v's"))
  print(mat, quote=F)
  cat('\n Best sumabs value : ', x$bestsumabsu, "\n")
  cat('\n Lambda Used : ', x$lambda, "\n")
  if(x$upos) cat("Elements of u constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("Elements of u constrained to be negative.", fill=TRUE)
}
