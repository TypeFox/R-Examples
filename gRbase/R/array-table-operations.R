tablePerm <- function(x, perm, resize=TRUE, keep.class=FALSE){
  # Like aperm() but perm can be dimnames
  if (missing( perm )){
    perm <- integer(0)
    return(aperm.default( x, perm, resize ))
  }

  if (is.character( perm )){
    perm <- match(perm, names(dimnames( x )))
    if ( any( is.na( perm )))
      stop("Invalid permutation...")
  }
  ans <- aperm.default( x, perm, resize )
  if (keep.class){
      class( ans ) <- oldClass( x )
  }
  ans
}

tableMult <- function(t1, t2){
  tableOp(t1,t2, op="*")
}

tableDiv <- function(t1, t2){
  tableOp(t1,t2, op="/")
}

tableOp <- function(t1, t2, op="*"){

  if (!is.array(t1)) {stop("'t1' is not an array")}
  if (!is.array(t2)) {stop("'t2' is not an array")}

  di1 <- dim(t1)
  di2 <- dim(t2)
  dn1 <- dimnames(t1)
  dn2 <- dimnames(t2)
  vn1 <- names(dn1)
  vn2 <- names(dn2)

  idx <- match(vn2, vn1)    ## location of variables in vn2 in vn1:
  idx.na <- is.na(idx)      ## logical of variables in {vn2\vn1}

  if (any(idx.na)){         ## If there are variables in {vn2 \ vn1}
    aug.vn <- vn2[idx.na]   ## Find those variables
    aug.di <- di2[idx.na]   ## - and their levels
    aug.dn <- dn2[idx.na]   ## - and their dimnames

    ## Create "augmented" table defined over (vn1, vn2\vn1) by repeating t1.
    pot1      <- rep.int(as.numeric(t1), prod(aug.di))
    vn.new    <- c(vn1, aug.vn)
    di.new    <- c(di1, aug.di)
    dn.new    <- c(dn1, aug.dn)
    dim(pot1)      <- di.new
    dimnames(pot1) <- dn.new
  } else {
    pot1   <- t1
    vn.new <- vn1
    di.new <- di1
    dn.new <- dn1
  }

  ## Find indices of vn2 in augmented table (vn1, vn2\vn1)
  vn2.idx    <- match(vn2, vn.new)
  ## Create perumation indices; first variables in vn2; then vn1\vn2
  perm  <-  c(vn2.idx, (1:length(vn.new))[-vn2.idx])

  if (op == "*") {
    pot1 <- as.numeric(aperm.default(pot1, perm, TRUE)) * as.numeric(t2)
  }
  else {
    pot1 <- as.numeric(aperm.default(pot1, perm, TRUE)) / as.numeric(t2)
    pot1[!is.finite(pot1)] <- 0
  }
  dim(pot1)      <- di.new[perm]
  dimnames(pot1) <- dn.new[perm]

  class(pot1) <- c("parray","array")
  pot1
}


.tableOp <- function(t1, t2, op="*"){

  if (!is.array(t1)) {stop("'t1' is not an array")}
  if (!is.array(t2)) {stop("'t2' is not an array")}

  op <- switch(op,
               "*"={`*`},
               "/"={`/`},
               "+"={`+`},
               "-"={`-`})

  di1 <- dim(t1)
  di2 <- dim(t2)
  dn1 <- dimnames(t1)
  dn2 <- dimnames(t2)
  vn1 <- names(dn1)
  vn2 <- names(dn2)

  idx <- match(vn2, vn1)   ## location of variables in vn2 in vn1:
  idx.na <- is.na(idx)      ## logical of variables in {vn2\vn1}

  if (any(idx.na)){         ## If there are variables in {vn2 \ vn1}
    aug.vn <- vn2[idx.na]   ## Find those variables
    aug.di <- di2[idx.na]   ## - and their levels
    aug.dn <- dn2[idx.na]   ## - and their dimnames

    ## Create "augmented" table defined over (vn1, vn2\vn1) by repeating t1.
    vn.new    <- c(vn1, aug.vn)
    di.new    <- c(di1, aug.di)
    dn.new    <- c(dn1, aug.dn)
    t1           <- rep.int(as.numeric(t1), prod(aug.di))
    dim(t1)      <- di.new
    dimnames(t1) <- dn.new
  } else {
    vn.new <- vn1
    di.new <- di1
    dn.new <- dn1
  }

  ## indices of vn2 in augmented table (vn1, vn2\vn1)
  vn2.idx    <- match(vn2, vn.new)
  ## Create perumation indices; first variables in vn2; then vn1\vn2
  perm  <-  c(vn2.idx, (1:length(vn.new))[-vn2.idx])

  tt1 <- op(aperm.default(t1, perm, TRUE), as.numeric(t2))
  if (identical(op, `/`))
    tt1[!is.finite(tt1)] <- 0
  dim(tt1)      <- di.new[perm]
  dimnames(tt1) <- dn.new[perm]
  class(tt1) <- c("parray","array")
  tt1
}

tableOp2 <- .tableOp2 <- function (t1, t2, op = `*`, restore = FALSE)
{
  if (!is.array(t1)){
    str( t1 )
    stop("'t1' is not an array")
  }
  if (!is.array(t2)){
    str( t2 )
    stop("'t2' is not an array")
  }

  vn1  <- names(dimnames(t1))
  vn2  <- names(dimnames(t2))

  ## indices of vn2 in vn1:
  vn2.idx   <- match(vn2, vn1)
  ## Create perumation indices; first variables in vn2; then vn1\vn2
  perm <- c(vn2.idx, (1:length(vn1))[-vn2.idx])

  pot1 <-
    if (restore) {
      zz    <- op(aperm.default(t1, perm, TRUE), as.numeric(t2))
      newvn <- c(vn2, vn1[-vn2.idx])
      perm2 <- match(vn1, newvn)
      aperm.default(zz, perm2, TRUE)
    } else {
      op(aperm.default(t1, perm, TRUE), as.numeric(t2))
    }
  if (identical(op, `/`))
    pot1[!is.finite(pot1)] <- 0
  pot1
}


tableSlice <-  function (x, margin, level, impose)
{
    if (is.null(margin))
        return(x)
    if (is.null(dimnames(x)))
        stop("'tableSlice' requires a structure with a dimnames attribute")

    dn    <- dimnames(x)
    vn    <- names(dn)

    if (is.character(margin)){
        mar.idx <- match(margin, vn)
        if (any((z<-is.na(mar.idx))))
            stop("Variable(s): ", margin[z], " do not exist in table...")
    } else {
        mar.idx <- margin
    }

    if (is.character(level)){
        lev.idx  <- rep(NA, length(level))
        for (kk in seq_along(margin)){
            lev.idx[kk] <- match(level[kk], dn[[mar.idx[kk]]])
        }
        if (any((z<-is.na(lev.idx))))
            stop("Level: ", level[z], " do not exist in table...")
    } else {
        lev.idx <- level
    }

    idx          <- vector("list", length(dim(x)))
    idx[]        <- TRUE
    idx[mar.idx] <- lev.idx
    ans <-do.call("[", c(list(x), idx))

    if (!missing(impose) && is.numeric(impose)){
        ans[] <- impose
    }
    ans <- array(ans, dim=sapply(dn[-mar.idx], length), dimnames=dn[-mar.idx])
    class(ans) <- c("parray","array")
    ans
}


## tableSlicePrim: Works only with margin and level being indices
tableSlicePrim <- function(x, mar.idx, lev.idx){
  idx         <- vector("list", length(dim(x)))
  idx[]       <-TRUE
  idx[mar.idx] <- lev.idx
  do.call("[", c(list(x), idx))
}

tableMargin <- function (x, margin, keep.class = FALSE)
{
##   cat("===== tableMargin =====\n")
##   print(as.data.frame.table(x));   print(margin)

    if (!is.array(x))
        stop("'x' is not an array")

    at <- attributes( x )
    di <- at[['dim']]
    dn <- at[['dimnames']]
    vn <- names( dn )

    if (length(margin)) {
        if( class(margin)=="formula" ){
            margin <- unlist(rhsf2list( margin ), use.names=FALSE)
        }
        if (is.character(margin)) {
            margin <- unique( margin )
            marg.idx <- match(margin, vn)
            if (any(is.na(marg.idx)))
                stop(sprintf("Variable(s): %s not in table ...\n",
                        toString( margin[is.na(marg.idx)] )) )
        }
        else {
            marg.idx <- margin
        }
        rest.idx <- (seq_along(vn))[-marg.idx]
        nr <- prod( di[marg.idx] )
        nc <- prod( di[rest.idx] )

        z <- rowSumsPrim(
            matrix(
                aperm.default(x, c(rest.idx, marg.idx), TRUE),
                nrow=nr, ncol=nc, byrow=TRUE))
        attributes(z) <- list(dim=di[marg.idx], dimnames=dn[marg.idx])
    } else {
        z <- sum(x)
    }
    if (keep.class)
        class(z) <- oldClass( x )
    return(z)
}


tableGetSliceIndex <- function(x, margin, level, complement=FALSE){
    di <- dim(x)
    dn <- dimnames(x)
    vn <- names(dn)
    nidx <- match(margin, vn)

    if ( any((z<-is.na(nidx))) ){
        stop(sprintf("margin %s not in table\n",
                     toString(margin[z])))
    }

    sidx <- unlist(lapply(seq_along(nidx),
                          function(i) {match(level[i], dn[[ nidx[i] ]])}),
                   use.names = FALSE)

    if (any((z<-is.na(sidx)))){
        stop(sprintf("level %s not in table\n", toString(level[z])))
    }

    out <- slice2entry(sidx, nidx, di)
    if (complement){
        (1:prod(di))[-out]
    } else {
        out
    }
}

tableSetSliceValue <- function(x, margin, level, complement=FALSE, value=0){
    idx <- tableGetSliceIndex(x, margin=margin, level=level, complement=complement)
    x[ idx ] <- value
    x
}



##
## Spring 2015
##


tabSlice2Entries <- function(x, slice, complement=FALSE){
    tabSlice2Entries_(x, names(slice), unlist(slice, use.names=FALSE), complement)
}

tabSlice2Entries_<- function(x, margin, level, complement=FALSE){
    di <- dim(x)
    dn <- dimnames(x)
    vn <- names(dn)

    marg <- match(margin, vn)

    lev <- unlist(lapply(seq_along(marg),
                         function(i) {match(level[i], dn[[ marg[i] ]])}),
                  use.names = FALSE)

    if ( any((z<-is.na(marg))) ){
        stop(sprintf("margin %s not in table\n", toString(margin[z])))
    }

    if (any((z<-is.na(lev)))){
        stop(sprintf("level %s not in table\n", toString(level[z])))
    }

    .tabSlice2Entries_(marg, lev, di, complement)
}

.tabSlice2Entries_ <- function(margin, level, dim, complement=FALSE){
    out <- slice2entry(level, margin, dim)
    if (complement){
        (1:prod(dim))[-out]
    } else {
        out
    }
}

## Multiplies entries in slice by value1 and all other entries by value2; if
## value1 or value2 are are NULL, then nothing happens
tabSliceMult <- function(x, slice, val=1, comp=0){
    if ( !is.null(val) ){
        idx <- tabSlice2Entries(x, slice)
        x[idx] <- x[idx] * val
    }
    if ( !is.null(comp) ){
        idx <- tabSlice2Entries(x, slice, complement=TRUE)
        x[idx] <- x[idx] * comp
    }
    x
}


tabCondProb <- function(tab, cond=NULL){
    if (!is.array( tab ))
        stop("tab must be an array")
    if (is.null(cond)){
        tab / sum( tab )
    } else {
        out <- tabDiv(tab, tabMarg(tab, cond))
        tabPerm(out, names(dimnames(tab)))
    }
}














