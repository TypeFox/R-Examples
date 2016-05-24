mmod_dimension <- function(object){
  .mmod_dimension(object$modelinfo$dlq, object$datainfo)
}

.mmod_dimension <- function(dlq, datainfo){
  .mmod_dimensionPrimitive(dlq, datainfo$cont.names, datainfo$disc.names, datainfo$disc.levels)
}

.mmod_dimensionPrimitive <- function(dlq, cont.names, disc.names, disc.levels ){

  lll <- dlq$linear
  ddd <- dlq$discrete
  qqq <- dlq$quadratic

##   cat(".mmod_dimensionPrimitive\n")
##   cat(sprintf("CHK: linear generators: %s\n", toString(lll)))
##   cat(sprintf("CHK: cont.names: %s\n", toString(cont.names)))
   
  n.cont   <- length(cont.names)
  d.levels <- c(lapply(disc.levels, length), recursive=TRUE)
  
  ## discrete dimension
  ## ------------------
  ddd.num  <- lapply(ddd, match, disc.names)
  disc.dim <- .loglinDim(ddd.num, d.levels)
##   cat(sprintf("CHK: disc.dim=%i\n", disc.dim))
  
  ## linear dimension
  ## ----------------
  if (length(lll)==0){
    lin.dim <- n.cont
  } else {
    lin.dim <- 0
    for (ii in seq_along(cont.names)){
      cvar <- cont.names[ii]
      ##cat(sprintf("cvar=%s\n", cvar))
      ## Find those linear generators containing 'cvar'
      bbb <- unlist(lapply(lll, function(xx) any(match(xx, cvar))))
      ##print(bbb)
      if (is.na(any(bbb))){ ## If there are no such generators
        zzz <- 1
      } else {
        ## Find discrete generators corresponding to 'cvar'
        ccc <- lapply(lll[which(bbb)], setdiff, cvar)        
        ## A linear generator may consist of a 'cvar' itself; such generator must be removed
        ccc <- ccc[ sapply(ccc,length) >0 ]
        ##print(ccc)
        if (length(ccc)>0){
          ccc.num <- lapply(ccc, match, disc.names)
          ##print(ccc.num)
          ##print(d.levels)
          zzz     <- .loglinDim(ccc.num, d.levels) + 1
          ##print(zzz)
        }
      }
      ##cat(sprintf("CHK:  .. cvar=%10s dim=%i\n", cvar, zzz))
      lin.dim <- lin.dim + zzz    
    }
  }
  ##cat(sprintf("CHK: lin.dim=%i\n", lin.dim))

  ## quadratic dimension
  ## -------------------
  uuu      <- glist2adjMAT(qqq)
  quad.dim <- sum(uuu[upper.tri(uuu)]) + nrow(uuu)
  ##cat(sprintf("CHK: quad.dim=%i\n", quad.dim))

  ## total dimension
  ## ---------------
  mod.dim <- disc.dim + lin.dim + quad.dim
  
  ## saturated model dimension
  ## -------------------------
  sat.dim <- prod(d.levels)-1 +  ## discrete 
    n.cont * prod(d.levels) +    ## linear 
      n.cont*(n.cont+1)/2        ## quadratic
  
  ## independence model dimension
  ## ----------------------------
  i.dim <- sum(disc.dim - 1) + 2*n.cont
  
  c(mod.dim=mod.dim, sat.dim=sat.dim, i.dim=i.dim, df=sat.dim-mod.dim, idf=mod.dim-i.dim)
}


##       browser()
##       idx  <- any(unlist(lapply(lll, function(xx) any(match(xx, cvar)))))
##       cat(sprintf("  .. cvar=%s\n", cvar))
##       if (is.na(idx)){
##         zzz <- 1
##       } else {
##         aaa  <- lll[idx][!unlist(lapply(lll[idx], is.null))]      
##         print(aaa)
##         aaa  <- lapply(aaa, function(xx) xx[-1])
##         print(aaa)
##         if (length(aaa)>0){
##           aaa.num <- lapply(aaa, match, disc.names)
##           print(aaa.num)
##           zzz     <- .loglinDim(aaa.num, disc.dim) + 1
##         }
##       }



.loglinDim <- function(glistNUM, dtab){

  #print(glistNUM); print (dtab)
  .subsets <- function(x) {
    y <- list(vector(mode(x), length = 0))
    for (i in seq_along(x)) {
      y <- c(y, lapply(y, "c", x[i]))
    }
    y[-1]
  }

  max.g.size <- max(unlistPrim(lapply(glistNUM, length)))
  #max.g.size <- 1000
  ## cat("max.g.size:", max.g.size, "\n")

  if (max.g.size < 10)
    {
      zz    <- unlist(lapply(glistNUM, .subsets), recursive=FALSE)
      unzz  <- unique.default(zz)
    }
  else 
    {
      unzz     <- .subsets(glistNUM[[1]])
      base.idx <- unlistPrim(lapply(unzz, function(terms) sum(2^(terms - 1)) ))  
      if (length(glistNUM)>1){
        for (ii in 2:length(glistNUM)){
          tmp      <- .subsets(glistNUM[[ii]])
          tmp.idx  <- unlistPrim(lapply(tmp, function(terms) sum(2^(terms - 1)) ))    
          unzz     <- c(unzz,tmp[!(tmp.idx %in% base.idx)])
          base.idx <- c(base.idx, tmp.idx[!(tmp.idx %in% base.idx)])
        }
      }
    }

  ans <- 0
  for (jj in 1:length(unzz))
    ans <- ans + prod(dtab[unzz[[jj]]]-1)

  ans
}




## MIdimension <- function(object){

##   lll <- object$dlq$linear
##   ddd <- object$dlq$discrete
##   qqq <- object$dlq$quadratic
  
##   cont.names   <- object$cont.names
##   disc.names   <- object$disc.names
##   disc.levels  <- object$disc.levels
  
##   d.dim <- c(lapply(disc.levels, length), recursive=TRUE)
  
##   ## discrete dimension
##   ddd.num <- lapply(ddd, match, disc.names)
##   disc.dim <- .loglinDim(ddd.num, d.dim)
##   ## linear dimension
##   lin.dim <- 0
##   for (ii in seq_along(cont.names)){
##     cont <- cont.names[ii]
##     idx <- unlist(lapply(lll, function(xx) any(match(xx, cont))))
##     aaa<- lll[idx][!unlist(lapply(lll[idx], is.null))]
##     aaa<- lapply(aaa, function(xx) xx[-1])
##     #print(aaa)
##     if (length(aaa)>0){
##       aaa.num <- lapply(aaa, match, disc.names)
##       #print(aaa.num)
##       zzz <- .loglinDim(aaa.num, d.dim)
##       lin.dim <- lin.dim + zzz
##     }
##   }
##   ## quadratic dimension
##   uuu <- ugListMAT(qqq)
##   quad.dim <- sum(uuu[upper.tri(uuu)]) + nrow(uuu)
  
##   mod.dim <- disc.dim + lin.dim + quad.dim
  
##   ## saturated dimension
##   sat.dim <- (length(cont.names)+1)*(prod(d.dim)-1) + length(cont.names)*(length(cont.names)+1)/2

##   ## independence dimension
##   i.dim <- sum(d.dim - 1) + 2*length(cont.names)
  
##   c(mod.dim=mod.dim, sat.dim=sat.dim, i.dim=i.dim, df=sat.dim-mod.dim, idf=mod.dim-i.dim)
## }
