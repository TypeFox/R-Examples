####
#### Re-implementation of simulate() function - quite fast...
#### Bristol, March 2008
####

simulate.grain <- function(object, nsim=1, seed=NULL, ...){

  if (!object$isCompiled){
    ##cat("Compiling (and propagating) model ...\n")
    object <- compile(object, propagate=TRUE)
  } else {
    if (!object$isPropagated){
      ## cat("Propagating model...\n")
      object <- propagate(object)
    }
  }

  plist  <- object$equipot
  cqlist <- object$rip$cli
  splist <- object$rip$sep

  ## Init
  ans           <- matrix(0, nrow=nsim, ncol=length(nodeNames(object)))
  colnames(ans) <- nodeNames(object)

  ctab  <- plist[[1]]
  res   <- simulateArray(x=ctab, nsim=nsim)
  ans[,colnames(res)] <- res

  ## Iterate
  if (length(cqlist)>1){
    for (ii in 2:length(cqlist)){
      ctab <- plist[[ii]]
      vn   <- names(dimnames(ctab))
      sp   <- splist[[ii]] ## What we condition on
      if (length(sp)>0){
        mtab <- tableMargin(ctab, sp)
        ctab <- tableOp2(ctab, mtab, `/`)
      }
      rr   <- setdiff(vn, sp) ## Variables to be simulated
      ##cat("r:", rr, "s:", sp, "\n")
      if (length(sp)){
        spidx <- match(sp, vn)
        res   <- matrix(0, nrow=nsim, ncol=length(rr))
        colnames(res) <- rr
        un    <- ans[, sp, drop=FALSE]
        ##cat("un:\n"); print(un)
        vals  <- unique(un)
        sc    <- cumprod(apply(vals, 2, max) )
        sc    <- c(1,sc)[1:length(sc)]
        key   <- ((un-1) %*% sc)+1
        ##cat(sprintf("key=%s\n", toString(key)))  #browser()
        for(kk in unique(key)){
          nn   <- sum(kk==key)
          idx  <- un[match(kk, key),]
          res[kk==key,] <- simulateArray(ctab, nsim=nn, margin=spidx, value.margin=idx)
        }
      } else {
        res <- simulateArray(x=ctab, nsim=nsim)
      }
      ans[,colnames(res)] <- res
    }
  }

  ns <- nodeStates(object)
  vn <- colnames(ans)
  aaa <- vector("list", ncol(ans))
  names(aaa) <- vn
  for (jj in 1:ncol(ans)){
    aaa[[jj]] <- factor(ans[,jj], levels=seq(ns[[jj]]))
    levels(aaa[[jj]]) <- ns[[jj]]
  }
  aaa <- as.data.frame(aaa)
  names(aaa) <- vn
  aaa
}



##   ans <- as.data.frame(ans)
##   vn <- names(ans)

##   for (jj in 1:ncol(ans)){
##     #match(vn[jj], names(ns))
##     ans[,jj] <- factor(ans[,jj], levels=seq(ns[[jj]]))
##     levels(ans[,jj]) <- ns[[jj]]
##   }

  #return(ans)


      ##cat(sprintf("vn=%s sp=%s\n", toString(vn), toString(sp)))
      ##cat("ctab:\n");  print(ctab)
      ##cat("mtab:\n"); print(mtab)
      ##cat("ctab (updated):\n"); print(ctab)
