#########################################################
##
## Find dimension of general log-linear model
##
## 'glistNUM' : generating class represented as numbers
## 'nlev'     : corresponding levels of the variables
##
#########################################################

loglinGenDim <- function(glist, tableinfo){

  if (is.table(tableinfo)){
    dtab <- c(lapply(dimnames(table), length), recursive=TRUE)
  } else {
    if (is.list(tableinfo)){
      dtab <- c(lapply(tableinfo, length), recursive=TRUE)
    } else {
      dtab <- tableinfo
    }
  }
  .loglinGenDim(glist, dtab)
}

.loglinGenDim <- function(glist, dtab){
  
  .subsets <- function(x) {
    y <- list(vector(mode(x), length = 0))
    for (ii in seq_along(x)) {
      y <- c(y, lapply(y, "c", x[ii]))
    }
    y[-1]
  }

  max.g.size <- max(unlistPrim(lapply(glist, length)))
  #cat("max.g.size:", max.g.size, "\n")

  if (max.g.size < 10)
    {
      #print(lapply(glist, .subsets))
      zz    <- unlist(lapply(glist, .subsets), recursive=FALSE)
      unzz  <- unique.default(zz)
      #print(unzz)
      
    }
  else 
    {
      unzz <- .subsets(glist[[1]])
      base.idx  <- unlistPrim(lapply(unzz, function(terms) sum(2^(terms - 1)) ))  
      if (length(glist)>1){
        for (ii in 2:length(glist)){
          tmp      <- .subsets(glist[[ii]])
          tmp.idx  <- unlistPrim(lapply(tmp, function(terms) sum(2^(terms - 1)) ))    
          unzz     <- c(unzz,tmp[!(tmp.idx %in% base.idx)])
          base.idx <- c(base.idx, tmp.idx[!(tmp.idx %in% base.idx)])
        }
      }
    }


  ans <- 0
  for (jj in 1:length(unzz)){
    inc <- prod(dtab[unzz[[jj]]]-1)
    ans <- ans + inc
  }
  ans
}


## Find dimension of decomposable model
## (with or without dimension adjustment for sparsity)
##
## 'cliq', 'seps' are cliques and separators (can be found from rip() function)
## 'table' can be either an array or a vector of levels with named components
##


loglinDecDim <- function(glist, tableinfo, adjust=TRUE){
  rr <- ripMAT(glist2adjMAT(glist))
  .loglinDecDim(rr$cliques, rr$separators, tableinfo=tableinfo, adjust=adjust)
}

.loglinDecDim <- function(cliq, seps, tableinfo, adjust=TRUE){

  if (!adjust){
    if (is.array(tableinfo))
      vlev <- c(lapply(dimnames(tableinfo), length), recursive=TRUE)
    else
      vlev <- tableinfo
    
    ## Without adjustment of dimension for sparsity
    npar <- prod(vlev[cliq[[1]]])-1
    if (length(cliq)>1){
      for (ii in 2:length(cliq)){
        dimC <- prod(vlev[cliq[[ii]]])-1
        dimS <- prod(vlev[seps[[ii]]])-1
        npar <- npar +  dimC - dimS  
      }
    }
  } else { 
    ## cat("With adjustment of dimension for sparsity\n")
    if (!is.array(tableinfo))
      stop("Model dimension adjusted for sparsity requires an array\n")
    tm1  <- tableMargin(tableinfo, cliq[[1]])
    npar <- sum(1*(tm1>0))-1    
    if (length(cliq)>1){
      for (ii in 2:length(cliq)){
        tm1  <- tableMargin(tableinfo, cliq[[ii]])
        tm2  <- tableMargin(tm1, seps[[ii]])
        dimC <- sum(1*(tm1>0))-1
        dimS <- sum(1*(tm2>0))-1
        npar <- npar + dimC - dimS  
      }
    }
  }
  return(npar)









## Extract from loglin().
## Just used for checking purposes
##
.getDimkh <- function(glist, dtab){

  subsets <- function(x) {
    y <- list(vector(mode(x), length = 0))
    for (i in seq_along(x)) {
      y <- c(y, lapply(y, "c", x[i]))
    }
    y[-1]
  }

  nvar <- length(dtab)
  df <- rep.int(0, 2^nvar)
  for (k in seq_along(glist)) {
    terms <- subsets(glist[[k]])
    for (j in seq_along(terms)){ 
      df[sum(2^(terms[[j]] - 1))] <- prod(dtab[terms[[j]]] -  1)
    }
   } 
  sum(df)
}





}




##   unzz <- subsets(glistNUM[[1]])
##   if (length(glistNUM)>1){
##     for (ii in 2:length(glistNUM)){
##       unzz <- unique(c(unzz, subsets(glistNUM[[ii]])))
##     }
##   }


## .getDim <- function(glistNUM, nlev, ind=0){
##   b <- glistNUM[[1]]
##   if (length(b)==0)
##     ans <- dimb <- 0 
##   else
##     ans <- dimb <- prod(nlev[b])-1
##   ##   if (length(b)>0){
##   ##     cat("glistNUM:\n"); 
##   ##     lapply(glistNUM, function(gg) cat(paste(paste(rep(" ",ind), collapse=''),paste(gg, collapse=" "),collapse=""),"\n"))
##   ##   }

##   if (length(glistNUM)>1){
##     B <- glistNUM[-1]
##     Bstar <- lapply(B, function(a) removeRedundant(intersectPrim(b, a)))
##     ans <- dimb + .getDim(B, nlev, ind+1) - .getDim(Bstar, nlev,ind+1)
##   }
##   return(ans)
## }

## ddd <- function(glistNUM, nlev, ind=0){
##   dd <- prod(nlev[glistNUM[[1]]])-1

##   if (length(glistNUM)>1){
##     for (ii in 2:length(glistNUM)){
##       bb <- glistNUM[[ii]]
##       B  <- glistNUM[1:(ii-1)]
##       Bstar <- lapply(B, function(a) removeRedundant(intersectPrim(bb, a)))
##       dd <- dd + prod(nlev[bb])-1 - ddd(Bstar, nlev)
##     }
##   }

##   dd
  

## }
