varNames.cptable <- function(x){
    x$vpa
}

valueLabels.cptable <- function(x){
    out <- list(x$levels)
    nam <- x$vpa
    names(out)<- x$vpa[1]
    out
}

varNames.table <- function(x){
    names(dimnames(x))
}

valueLabels.table <- function(x){
    dimnames(x)
}



cptable <- function(vpar, levels=NULL, values=NULL, normalize=TRUE,  smooth=0 ){
  vpa  <- c(.formula2char(vpar))
  ans  <- list(vpa=vpa, values=values, normalize=normalize, smooth=smooth, levels=levels)
  class(ans) <- "cptable"
  ans
}

ortable <- function(v, pa1=c(TRUE,FALSE), pa2=c(TRUE,FALSE), levels ){

  vpa <- c(.formula2char(v))

  if (length(vpa)!=3)
    stop("Must have exactly two parents!")
  lpa1 <- length(pa1)
  lpa2 <- length(pa2)
  z <- rep(pa1,lpa2) | rep(pa2,each=lpa1)
  pp <- array(c(z, !z),c(2,lpa2,lpa1))
  values <- as.numeric(aperm(pp, c(3,1,2)))
  ans <- list(vpa=vpa, values=values, normalize=FALSE, smooth=0, levels=levels)
  class(ans) <- "cptable"
  return(ans)
}

andtable <- function(v, pa1=c(TRUE,FALSE), pa2=c(TRUE,FALSE), levels ){

  vpa <- c(.formula2char(v))

  if (length(vpa)!=3)
    stop("Must have exactly two parents!")
  lpa1 <- length(pa1)
  lpa2 <- length(pa2)
  z    <- rep(pa1,lpa2) & rep(pa2,each=lpa1)
  pp   <- array(c(z, !z),c(2,lpa2,lpa1))
  values    <- as.numeric(aperm(pp, c(3,1,2)))
  ans       <- list(vpa=vpa, values=values, normalize=FALSE, smooth=0, levels=levels)
  class(ans) <- "cptable"
  return(ans)
}


print.cptable <- function(x,...){
  cat(sprintf("{v,pa(v)}      : %s\n", toString(x$vpa)))
  cat(sprintf("levels of v    : %s\n", toString(x$levels)))
  cat(sprintf("values         : %s\n", toString(x$values)))
  cat(sprintf("normalize=%s, smooth=%f\n", x$normalize, x$smooth))
  return(invisible(x))
}


print.cptable <- function(x,...){
    v <- x$values
    dim(v) <- c(length(x$levels),length(v)/length(x$levels))
    rownames(v) <- x$levels
    colnames(v) <- rep(NA, ncol(v))
    cat(sprintf("{v,pa(v)} :"))
    str(x$vpa)
    print(v)
    ## cat(sprintf("{v,pa(v)}      : %s\n", toString(x$vpa)))
  ## cat(sprintf("levels of v    : %s\n", toString(x$levels)))
  ## cat(sprintf("values         : %s\n", toString(x$values)))
  ## cat(sprintf("normalize=%s, smooth=%f\n", x$normalize, x$smooth))
  return(invisible(x))
}

