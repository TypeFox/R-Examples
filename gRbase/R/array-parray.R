######################################################
####
#### A parray is an array in the sense that it is a
#### vector with a dim and a dimnames attribute.
####
######################################################

parray <- function(varNames, levels, values=1, normalize="none", smooth=0){

  normalize <- match.arg(normalize, choices=c("none","first","all"))
  varNames  <- rhsFormula2list(varNames)[[1]]
  if (smooth>0){
    values <- values + smooth
  }

  dn   <- makeDimNames(varNames, levels)
  nlev <- unlist(lapply(dn, length))
  ans  <- array(values, dim=nlev, dimnames=dn)

  ## Normalize if requested
  switch(normalize,
         "first" = {
           ##cat("first\n")
           if (length(nlev)>1){
             tmp   <- matrix(ans, ncol=dim(ans)[1], byrow=TRUE)
             ans[] <- t.default(tmp/rowSumsPrim(tmp))
           } else {
             ans <- ans / sum(ans)
           }},
         "all"  = {
           ans <- ans / sum(ans)
         },
         "none" = { } )
  class(ans) <- c("parray","array")
  return(ans)
}

as.parray  <- function(values, normalize="none", smooth=0){

  normalize <- match.arg(normalize, choices=c("none","first","all"))

  if (!inherits(values, c("array","matrix","integer","double","table"))){
    stop("arg must be array, matrix, table, integer or double\n")
  }

  if (smooth>0){
    values <- values + smooth
  }

  if (is.null(dimnames(values))){
    if (!is.null(dim(values)))
      nLevels <- dim(values)
    else
      nLevels <- length(values)

    varNames <- paste("V", 1:length(nLevels),sep='')
    dimnames <- makeDimNames(varNames, nLevels)
    ans <- array(values, dim = nLevels, dimnames = dimnames)
  } else {
    ans <- values
  }
  class(ans) <- c("parray","array")

  switch(normalize,
    "first"={
      if (length(dim(ans))>1){
        marg  <- 2:length(dim(ans))
        ma    <- apply(ans, marg, sum)
        ans   <- sweep(ans, marg, ma, "/")
      } else {
        ans <- ans / sum(ans)
      }
    },
    "all"={ans <- ans / sum(ans)
    },
    "none"={}
    )
  attr(ans, "call") <- NULL
  return(ans)
}

data2parray <- function(data, varNames=NULL, normalize="none", smooth=0){
  cls <- match(class(data), c("data.frame","table", "xtabs", "matrix"))[1]
  if (is.na(cls)){
    stop("'data' must be one of  dataframe, table, xtabs, matrix")
  }

  .set.varNames <- function(varNames, dataNames){
    if (is.null(varNames)){
      if (is.null(dataNames))
        stop("'data' has no variable names")
      varNames <- dataNames
    } else {
      if (class(varNames) %in% c("formula", "character")){
        varNames <- rhsf2list(varNames)[[1]]
      }
    }
    varNames
  }

  switch(as.character(cls),
         "1"={
           dataNames <- names(data)
           varNames <- .set.varNames(varNames, dataNames)
           val  <- xtabs(~., data = data[, varNames, drop = FALSE])
         },
         "2"=, "3"=, "4"={
           dataNames <- names(dimnames(data))
           varNames <- .set.varNames(varNames, dataNames)
           val  <- tableMargin(data, varNames)
         }
         )
  attr(val, "call") <- NULL
  res <- as.parray(val, normalize = normalize, smooth = smooth)
  res
}



makeDimNames <- function(varNames, levels, sep=''){
  if (missing(varNames) || is.null(varNames))
    return(lapply(levels, seq))
  mapply(function(vv,ll){
    if (!is.character(ll)){
      if (length(ll)==1){
        ll <- 1:ll
      }
      ll <- paste(vv,ll,sep="")
    }
    ll
  }, varNames, levels, SIMPLIFY=FALSE)
}

makeDimNames <-
function(varNames, levels, sep=''){
    if (missing(varNames) || is.null(varNames))
        return(lapply(levels, seq))
    if (length(varNames) != length(levels))
        stop("'varNames' and 'levels' must have the same length")
    if (is.list(levels)){
        names(levels) <- varNames
        return( levels )
    }
    out <- lapply(seq_along(varNames), function(i)
                  {
                      ll <- levels[[ i ]]
                      vv <- varNames[ i ]
                      if (!is.character(ll)){
                          if (length(ll)==1){
                              ll <- 1:ll
                          }
                          ll <- paste(vv,ll,sep="")
                      }
                      ll
                  })
    names(out) <- varNames
    out
}

varNames.array     <- function(x) names(attr(x,"dimnames"))
nLevels.array      <- function(x) dim(x)
valueLabels.array  <- function(x) attr(x,"dimnames")

varNames.parray    <- function(x) names(attr(x,"dimnames"))
nLevels.parray     <- function(x) dim(x)
valueLabels.parray <- function(x) attr(x,"dimnames")

print.parray  <- function(x,...){
  class(x)<-NULL
  print(x)
  invisible(x)
}




















## ptable <- function(varNames, levels, values=1, normalize=c("none","first","all"), smooth=0){

##   normalize <- match.arg(normalize, choices=c("none","first","all"))
##   varNames  <- rhsFormula2list(varNames)[[1]]

##   if (is.list(levels)){
##     dimnames        <- levels
##     names(dimnames) <- varNames
##     levels          <- sapply(dimnames, length)
##   } else {
##     dimnames <- makeDimNames(varNames, levels)
##   }

##   if (smooth>0){
##     values <- values + smooth
##   }

##   ans <- array(values, dim=levels, dimnames=dimnames)

##   ## Normalize if requested
##   switch(normalize,
##          "first" = {
##            if (length(dim(ans))>1){
##              marg  <- 2:length(dim(ans))
##              ma2   <- tableMargin(ans, marg)
##              ans   <- tablePerm(.tableOp2(ans, ma2, op=`/`), names(dimnames(ans)))
##            } else {
##              ans <- ans / sum(ans)
##            }},
##          "all"  = {ans <- ans / sum(ans) },
##          "none" = { } )
##   class(ans) <- "ptable"
##   return(ans)
## }

## ## Create list with dimension names
## ##
## makeDimNames <- function(varNames, levels, sep=''){
##   if (missing(varNames) || is.null(varNames))
##     return(lapply(levels, seq))
##   lev <- lapply(levels, function(a) c(1:a))
##   mapply(function(n,v) paste(n,v,sep=sep), varNames, lev, SIMPLIFY=FALSE)
## }
## Accessors
##












