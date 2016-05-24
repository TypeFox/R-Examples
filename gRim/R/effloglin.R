##
## Implementation of efficient IPS algorith for log-linear models
## Based on fitting model using a grain structure
##
## Argument names are chosen so as to match those of loglin()
##

effloglin <- function(table, margin, fit=FALSE, eps=0.01, iter=20, print=TRUE){
  
  amat <- ugList(margin, result="matrix")
  vn   <- colnames(amat)
  tri  <- triangulateMAT(amat)
  rip  <- ripMAT(tri)

  cliq     <- rip$cliques
  len.cliq <- length(cliq)

  ## "host clique" for each generator
  ##
  ghost   <- rep(NA, length(margin))
  seqcliq <- seq_along(cliq)
  for (kk in 1:length(margin)){
    ##cat("kk:", kk,"\n")
    gg <- margin[[kk]]
    for (ii in seqcliq){
      ##cat ("ii", ii, "\n")
      zz <- match(gg, cliq[[ii]])
      if (!any(is.na(zz))){
        ghost[kk] <- ii
        break
      } 
    }
  }
    
  if (is.array(table)){
    Nobs   <- sum(table)
    stlist <- lapply(margin, function(xx) {tableMargin(table, xx)})    
  } else {
    Nobs   <- sum(table[[1]])
    stlist <- table
  }

  zzz       <- unlist(lapply(stlist, dimnames), recursive=FALSE)
  vl        <- zzz[uniquePrim(names(zzz))]
  pot.list  <- lapply(cliq, function(cq)
                      parray(cq, levels=vl[cq], values=1, normalize="all"))
##   cat("effloglin\n")
##   print(as.data.frame.table(pot.list[[1]]))
  
  ## ## Potential list over cliques
  ## Clique marginals
  prob.list  <- propagateLS(pot.list, rip, initialize=TRUE)        

  itcount  <- 1L
  logL     <- 0
  zzz <- vector("numeric", length(margin))
  repeat{
    cat(sprintf("---------- iteration: %i -----------\n", itcount))
    for (ss in seq_along(margin)){
      gg      <- margin[[ss]]
      st      <- stlist[[ss]]
      cq      <- cliq[[ghost[ss]]]
      cq.idx  <- ghost[ss]      
      cpot    <- prob.list[[cq.idx]]
      ##adjust  <- tableOp(st, tableMargin(cpot, gg)*Nobs, "/")
      
      tm      <- tableMargin(cpot, gg)*Nobs
      adjust  <- st / tm
      zzz[ss] <- max(abs(log(adjust)))
      ##zzz[ss] <- max(abs(st-tm))
      logL    <- logL + sum(st * log(adjust))
      ##pot.list[[cq.idx]] <- tableOp(pot.list[[cq.idx]], adjust, "*")
      pot.list[[cq.idx]] <- tableOp2(pot.list[[cq.idx]], adjust, `*`)
      prob.list          <- propagateLS(pot.list, rip, initialize=TRUE)
    }
    
    if (print)
      cat("max deviation (obs-fitted):", max(zzz), "\n")
    if (max(zzz)<eps || itcount>=iter)
      break()
    itcount <- itcount + 1L
  }
   
  vl    <- unlist(lapply(stlist, dimnames), recursive=FALSE)[vn]
  nlev  <- unlistPrim(lapply(vl, length))  
  gn    <- lapply(margin, match, vn)
  nparm <- .loglinGenDim(gn, nlev)
  df    <- prod(nlev) - 1 - nparm
  
  ans <- list(potlist=pot.list, margin=margin, vn=vn, rip=rip, ghost=ghost,
              stlist=stlist, logL=logL, nparm=nparm, df=df)
  
### Create full joint:
  if (fit){
    pjoint <- prob.list[[1]]
    if (length(prob.list)>1){
      for (ii in 2:length(prob.list)){
        pjoint <- tableOp(pjoint, tableOp(prob.list[[ii]],
                                          tableMargin(prob.list[[ii]], rip$sep[[ii]]),
                                          "/"),"*")
      }
    }
    pjoint <- tablePerm(pjoint, vn)*Nobs    
    ans <- c(ans, list(fit=pjoint))
  }  
  ## class(ans) <- "effloglin"
  return(ans)
}





