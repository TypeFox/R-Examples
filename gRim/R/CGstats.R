#####################################################################
###
### Statistics for conditional Gaussian distribution
###
#####################################################################

CGstats <- function(object, varnames=NULL, homogeneous=TRUE, simplify=TRUE){
  UseMethod("CGstats")
}


CGstats.data.frame <- function(object, varnames=NULL, homogeneous=TRUE, simplify=TRUE){

  if (is.null(varnames)){
    varnames <- names(object)
  }
  
  if (is.numeric(varnames)){
    use.idx <- varnames
  } else {
    use.idx <- match(varnames, names(object))
  }
  
  zzz <- unlist(lapply(object, is.numeric))

  cont.idx <- intersect(which(zzz),  use.idx)
  disc.idx <- intersect(which(!zzz), use.idx)
  
  cont.names <- names(zzz)[cont.idx]
  disc.names <- names(zzz)[disc.idx]

  CGstats_internal(object, disc.names, cont.names, homogeneous, simplify)
}

CGstats_internal <- function(object, disc.names=NULL, cont.names=NULL, homogeneous=TRUE, simplify=TRUE){

  if (length(cont.names)==0)
    {
      xt    <- xtabs(~., data=object[,disc.names,drop=FALSE])
      ans   <- list(n.obs=xt)
      disc.levels <- dim(xt)
    }
  else
    {
      if (length(disc.names)==0)
        {
          ans <- cov.wt(object[,cont.names,drop=FALSE], method="ML")
          disc.levels <- NULL
        }
      else
        {
          ans <- .cov.wt.by(disc.names, cont.names, object)
          disc.levels <- dim(ans$n.obs)
          
          if (homogeneous) {
            pp   <- ans$n.obs / sum(ans$n.obs)
            CC   <- ans$cov[[1]]
            CC[] <- 0
            for (ii in 1:length(ans$cov)){
              CC <- CC + ans$cov[[ii]]*pp[ii]
            } 
            ans$cov <- CC 
          }
          
          if (simplify){
            ans$center <- t(do.call(rbind, ans$center))
            rownames(ans$center) <- cont.names
            if (!homogeneous){
              ans$cov <- t(do.call(rbind, lapply(ans$cov, as.numeric)))
            }
          }
        }
    }
  res  <- c(ans, list(cont.names=cont.names, disc.names=disc.names, disc.levels=disc.levels))
  ##class(res) <- "CGstats"
  return(res)
 }

print.CGstats <- function(x,...){
  print.default(x[1:3])
  return(invisible(x))
}


.cov.wt.by <- function(disc.names, cont.names, data){
  
  sss   <- .split.by(disc.names, indata=data)
  
  SS.mean  <- SS.cov   <- vector("list", length(sss$uniqval))
  names(SS.mean)  <- names(SS.cov) <- sss$uniqval
  
  for (ii in 1:length(sss$uniqval)){
    xx <- data[sss$uniqval[ii]==sss$facstr,cont.names,drop=FALSE]
    zz <- .cov.wt(xx, method="ML" )
    SS.mean[[ii]] <- zz$center
    SS.cov[[ii]]  <- zz$cov
  }
  
  ans <- list(n.obs=sss$xt, center=SS.mean, cov=SS.cov)
  return(ans)

}

.cov.wt <- function(xx,method="ML"){
  xx <- as.matrix(xx)
  n.obs <- nrow(xx)
  center<-colSumsPrim(xx) / n.obs
  zz  <- xx - rep(center,each=nrow(xx))
  if (method=="ML")
    NNN <- n.obs
  else
    NNN <- n.obs-1

  ccc <- crossprod(zz) / NNN
  ans <- list(cov=ccc, center=center, n.obs=n.obs)
  ans
}

## Works only for heterogeneous, mixed variables and simplify=TRUE
.extendCGstats <- function(CGstats) {
  n.i    <- as.numeric(CGstats[['n.obs']])
  N      <- sum(n.i)
  Q      <- nrow(CGstats[['center']])
    
  SSD    <- matrix(rowSumsPrim(.colmult(n.i, CGstats[['cov']])), nrow=Q)
  S      <- SSD/N
    
  mmm   <- CGstats[['center']]
  ttt   <- .colmult(n.i, mmm)
  quad  <- ttt %*% t(mmm)
  SS    <- SSD + quad
  
  ans <- c(CGstats, list(N=N, SSD=SSD, SS=SS))
  ##class(ans) <- "CGstats"
  ans
}

.split.by <- function(disc.names, indata, drop=FALSE){

  xt   <- xtabs(~., indata[,disc.names,drop=FALSE])
  ooo  <- as.data.frame.table(xt, useNA="always")[,1:length(disc.names),drop=FALSE]

  ooostr <- .dfcols2namevec(ooo)
  facstr <- .dfcols2namevec(indata[,disc.names,drop=FALSE])

  ans <- list(facstr=facstr, uniqval=ooostr, xt=xt)
  return(ans)
  }


.dfcols2namevec <- function(ooo,sep='|'){
  ans <- as.character(ooo[,1])
  if (ncol(ooo)>1){
    for (jj in 2:ncol(ooo)){
      ans <- paste(ans, ooo[,jj],sep=sep)
    }
  }
  return(ans)
}


