setClass("paras",
         representation = representation(vals = "numeric", pnames="character")
         )

setClass("MB", ## "MB" == "multiplicative binomial"
         representation = representation(counts = "matrix", m="numeric", pnames="character")
         )

setGeneric("pnames",function(x){standardGeneric("pnames")})
setGeneric("pnames<-",function(x, value){standardGeneric("pnames<-")})

setMethod("pnames","paras",function(x){x@pnames})
setMethod("pnames","MB",function(x){x@pnames})

setGeneric("getVals",function(x){standardGeneric("getVals")})
setMethod("getVals","paras",function(x){x@vals})

setGeneric("counts",function(x){standardGeneric("counts")})
setMethod("counts","MB",function(x){
  out <- x@counts
  colnames(out) <- pnames(x)
  return(out)
} )

setGeneric("getM",function(x){standardGeneric("getM")})
setMethod("getM","MB",function(x){
  out <- x@m
  names(out) <- pnames(x)
  return(out)
} )

### There are no occurrences of "@" below this line

setOldClass("Oarray")

## paras first, then MB:

setMethod("pnames<-","paras",function(x,value){
  paras(getVals(x),pnames=value)
} )

setMethod("pnames<-","MB",function(x,value){
MB(counts(x),getM(x),pnames=value)
} )

setGeneric("p",function(x){standardGeneric("p")})
setMethod("p","paras",function(x){
  jj <- getVals(x)[seq_len(length(x)-1)]
  jj <- c(jj,1-sum(jj))
  pn <- pnames(x)
  if(length(pn)>0){
    names(jj) <- pn
  }
  return(jj)
} )

setGeneric("p<-",function(x,value){standardGeneric("p<-")})

setMethod("p<-","paras",function(x,value){
  jj <- getVals(x)
  jj[seq_len(length(x)-1)] <- value
  paras(jj, pnames=pnames(x))
} )

setGeneric("theta",function(x){standardGeneric("theta")})
setMethod("theta","paras",function(x){
  k <- length(x)
  jj <- matrix(1,k,k)
  jj[upper.tri(jj,FALSE)] <-  getVals(x)[-seq_len(k-1)] 
  
  rownames(jj) <- pnames(x)
  colnames(jj) <- pnames(x)
  jj
} )

setGeneric("theta<-",function(x,value){standardGeneric("theta<-")})
setMethod("theta<-","paras",function(x,value){
  jj <- getVals(x)
  if(length(value)==1){
    jj[-seq_len(length(x)-1)] <- value
  } else {
    jj[-seq_len(length(x)-1)] <- value[upper.tri(value,FALSE)]
  }
  paras(jj, pnames=pnames(x))
} )

".fun" <- function(n){-0.5 + 0.5*sqrt(9+8*n)}
setMethod("length","paras",function(x){ .fun(length(getVals(x)))})

".print_paras_worker" <- 
function(x, ...){ # This function does the work

  jj <- theta(x)
  if(FALSE){
    jj <- round(jj,unlist(options("digits")))
    storage.mode(jj) <- "character"
    jj[lower.tri(jj,diag=FALSE)] <- "-"
    diag(jj) <- 1
  } else {
    jj[lower.tri(jj,diag=TRUE)] <- NA
  }
  
  return(list(p=p(x),theta=noquote(jj)))
}
    
".print_paras" <- function(x, ...){
  jj <- .print_paras_worker(x, ...)
  print(jj)
  return(invisible(jj))
}

setMethod("show", "paras",
          function(object){.print_paras(object)}
          )

".paras_valid" <- function(object){
  jj <- .fun(length(getVals(object)))
  if(!(abs(jj - round(jj)) < 0.01)){
    return("length not of form (n^2+n)/2-1")
  } else if (
             (length(pnames(object)) != 0 ) &
             (length(object) != length(pnames(object)) ) ){
    return("pnames not correct length")
  } else {
    return(TRUE)
  }
}

setValidity("paras", .paras_valid)

"paras" <- function(x, p, theta, pnames=character(0)){
  if(is.null(pnames)){pnames <- character(0)}
  if(!missing(x)){
    if(length(x)==1){
      return(new("paras",vals=c(rep(1/x,x-1), rep(1,x*(x-1)/2)),pnames=pnames))
    } else {
      return(new("paras", vals = x, pnames=pnames))
    }
  } else { 
    return(Recall(c(p,theta[upper.tri(theta,FALSE)])))
  }
}

setMethod("[", "paras", function(x, i, j, drop=TRUE){
  if(missing(j)){
    return(p(x)[i])
  } else {
    return(theta(x)[i,j,drop=drop])
  }
} )

setReplaceMethod("[",signature(x="paras"), function(x,i,j,value){
  if(missing(j)){
    temp_p <- p(x)[-length(x)]
    temp_p[i] <- value
    p(x) <- temp_p
    return(x)
  } else {
    temp_theta <- theta(x)
    temp_theta[i,j] <- value
    theta(x) <- temp_theta
    return(x)
  }
} )


## done paras, now MB:

".MB_valid" <- function(object){
  m <- getM(object)
  jj <- c(counts(object),m)

  pp <- sweep(counts(object),2,getM(object))
  if(any(jj-round(jj)>0.1)){
    return("non-integer entries")
  } else if(any(pp>0)){ # sic
    return('negative counts not allowed')
  } else if (any(jj<0)){
    return("negative entries")
  } else if(ncol(counts(object)) != length(m)){
    return("'m' wrong length")
  } else if(length(pnames(object)) != length(m)){
    return('pnames wrong length')
  } else {
    return(TRUE)
  }
}

setValidity("MB", .MB_valid)

"MB" <- function(dep, m, pnames=character(0)){

  if(inherits(dep,"Oarray")){
    return(.Oarray_to_MB(dep))
  }

  if(inherits(dep,"gunter_MB")){
    stop("problem in MB: method not implemented (yet)")
  }
  
  if(length(pnames)==0){
    if(is.null(colnames(dep))){
      if(!is.null(names(m))){
        pnames <- names(m)
      }
    } else {
      pnames <- colnames(dep)
    }
  }

  if(length(pnames)==0){pnames <- paste("Var", seq_along(m),sep="")}
  colnames(dep) <- NULL
  new("MB",counts=dep, m=m, pnames=pnames)
}

## getM() defined at the top of MM.R, using the 'at' symbol'

".print_MB" <- function(x, separator='   ', ...){
  jj <- counts(x)
  k <- length(getM(x))
  wanted <- c(t(cbind(matrix(seq_len(2*k),k,2,byrow=FALSE),2*k+1)))
  jjbar <- -sweep(jj,2,getM(x))
  colnames(jjbar) <- paste("n",colnames(jjbar),sep='')
  jj <- cbind(cbind(jj,jjbar," "=separator)[,wanted])
  print(noquote(jj))
  return(invisible(jj))
}

setMethod("show","MB",
          function(object){.print_MB(object)}
          )

"lmultinomial" <- function(x){lfactorial(sum(x))-sum(lfactorial(x))} 
 "multinomial" <- function(x){exp(lmultinomial(x))}
 
"MM_single" <- function(yrow,paras,givelog=FALSE){ # no normalizing constant; 'yrow' means a single row  of y.
  stopifnot(is.vector(yrow))
  M <- theta(paras)  # 'M' for 'Matrix', avoid theta nameclash
  p <- p(paras)
  M[lower.tri(M,TRUE)] <- 1 
  # multinomial(yrow) * prod(p^yrow) * exp(quad.form(log(theta),yrow))
  out <- lmultinomial(yrow) + sum(log(p)*yrow) + quad.form(log(M),yrow)
  if(givelog){
    return(out)
  } else {
    return(exp(out))
  }
} 

"dMM" <- function(Y,paras){ # includes normalizing constant
  p <- p(paras)
  theta <- theta(paras)
  MM_single(Y,paras)/NormC(sum(Y),paras) 
} 

"NormC" <- function(Y,paras,log=FALSE) {
  p <- p(paras)
  theta <- theta(paras)
  p <- p/sum(p) 
  jj <-   sum(apply(compositions(Y,length(p)),2,MM_single,paras)) 
  if(log){ 
    return(log(jj)) 
  } else { 
    return(jj) 
  } 
} 

"MM" <- function(y,n=NULL, paras){ # each row of 'y' is an observation
  p <- p(paras)
  theta <- theta(paras)
        
  if(minmax(rowSums(y))){ 
    return(MM_allsamesum(y,n,paras))
  } else { 
    return(MM_differsums(y,n,paras)) 
  } 
} 
 
"MM_allsamesum" <- function(y, n=NULL, paras){ # each row of 'y' is an observation; 'n' is the number of times that row was observed

  p <- p(paras)  # sum(p)=1
  theta <- theta(paras)
 
  theta[lower.tri(theta,TRUE)] <- 1 
  if(is.null(n)){n <- rep(1,nrow(y))} 
  Y <- sum(y[1,]) 
  N <- nrow(y) 
  f <- function(y){sum(lfactorial(y))}
  jj <-  NormC(Y,paras,log=TRUE)
      return( 
           -jj                           *   sum(n) 
           +lfactorial(Y)                *   sum(n) 
           -sum(apply(y,1,f)                   * n)
           +sum(rowSums(sweep(y,2,log(p),"*")) * n) 
           +sum(diag(quad.tform(log(theta),y)) * n)
           )
} 

"MM_allsamesum_A" <- function(y,paras){ # each row of 'y' is an observation
  p <- p(paras)
  theta <- theta(paras)
  stopifnot(minmax(rowSums(y))) 
  theta[lower.tri(theta,TRUE)] <- 1 
  Y <- sum(y[1,]) 
  N <- nrow(y)
  return(
         -N*NormC(Y,paras)
         +N*lfactorial(Y) 
         -sum(lfactorial(y)) 
         +sum(sweep(y,2,log(p),"*")) 
         +tr(quad.tform(log(theta),y)) 
         )
} 
 
"MM_support" <- function(paras, ss){ # 'ss' is a sufficient statistic;
                                     # use suffstats() to generate
                                     # this.  Observe that this is a
                                     # *support* function of paras.
  return(
         -NormC(ss$Y,paras,log=TRUE) * ss$nobs
         +sum(log(p(paras))          * ss$row_sums)
         +sum(log(theta(paras))      * ss$cross_prods)
         )
}

"MM_differsums" <- function(y, n=NULL, paras) { # each row of 'y' is an observation; 'n' is the number of times that row was observed
  p <- p(paras)
  theta <- theta(paras)
  
  if(is.null(n)){n <- rep(1,nrow(y))} 
  return(
         sum(n*( 
                -sapply(rowSums(y),NormC,paras,log=TRUE)
                +apply(y,1,lmultinomial)         
                +rowSums(sweep(y,2,log(p),"*"))   
                +diag(quad.tform(log(theta),y)) 
                ))
         )
} 

"MM_differsums_A" <- function(y, paras){  # here for completeness
  return(MM_differsums(y,n=NULL,paras))
}
 
"optimizer" <- function(y, n=NULL, start=NULL, method="nlm", printing=FALSE,  give_fit=FALSE,  ...){
  if(minmax(rowSums(y))){
    return(optimizer_allsamesum(y, n=n, start=start, method=method, printing=printing,  give_fit=give_fit,  ...))
  } else {
    return(optimizer_differsums(y, n=n, start=start, method=method, printing=printing,  give_fit=give_fit,  ...))
  }
}

"optimizer_allsamesum" <- function(y, n=NULL, start=NULL, method="nlm", printing=FALSE,  give_fit=FALSE,  ...){ 

  stopifnot(minmax(rowSums(y)))
  
  if(is.null(start)){ 
    start <- Lindsey(y,n)
  }

  jjn <- pnames(start)
  if(is(start,"paras")){
    start <- log(getVals(start))
  }

  ss <- suffstats(y,n)

  f <- function(xin){
    x <- exp(xin)
    jj <- paras(x)
    p <- p(jj)
    theta <- theta(jj)
    if(printing){
      print(paras(x,pnames=jjn))
    }
    if(any(p<0) | any(theta<0) ){ 
      out <- .Machine$double.xmax
    } else { 
      out <- -MM_support(paras(x),ss)
        ## out <- -MM_allsamesum(y,n=n,paras(x)) 
    }
    if(printing){
      print(noquote(paste("support is ",out)))
      print(noquote("-----"))
    }
    return(out) 
  } 
  if(method=="Nelder"){
    out <- optim(start, f, ...)
    short_out <- out$par
  } else if (method == "nlm") {
    out <- nlm(f,start,...)
    short_out <- out$estimate
  } else {
    print("method must be Nelder or nlm")
    stop()
  }

  jj <- paras(exp(short_out))
  pnames(jj) <- colnames(y)
  if(give_fit){
    return(list(fit=out, p=jj))
  } else {
    return(jj)
  }
  
}

"optimizer_differsums" <- function(y, n=NULL, start=NULL, method="nlm", printing=FALSE,  give_fit=FALSE,  ...){ 

  if(is.null(start)){ 
    start <- paras(ncol(y))
  }

  jjn <- pnames(start)
  if(is(start,"paras")){
    start <- log(getVals(start))
  }

  f <- function(xin){
    x <- exp(xin)
    jj <- paras(x)
    p <- p(jj)
    theta <- theta(jj)
    if(printing){
      print(paras(x,pnames=jjn))
    }
    if(any(p<0) | any(theta<0) ){ 
      out <- .Machine$double.xmax
    } else { 
      out <- -MM_differsums(y,n=n, paras(x))
    }
    if(printing){
      print(noquote(paste("support is ",out)))
      print(noquote("-----"))
    }
    return(out) 
  } 
  if(method=="Nelder"){
    out <- optim(start, f, ...)
    short_out <- out$par
  } else if (method == "nlm") {
    out <- nlm(f,start,...)
    short_out <- out$estimate
  } else {
    print("method must be Nelder or nlm")
    stop()
  }

  jj <- paras(exp(short_out))
  pnames(jj) <- colnames(y)
  if(give_fit){
    return(list(fit=out, p=jj))
  } else {
    return(jj)
  }  
} 

setGeneric("gunter",function(obs){standardGeneric("gunter")})

".gunter" <- function(obs){  ## slightly modified version of a function written by Bert Gunter, R-help, 16 October 2009.
  ## Try gunter(pollen) -- takes about a minute to run

  jj <- rowSums(obs)
  stopifnot(minmax(jj))
  tbl <- as.matrix(t(compositions(jj[1],ncol(obs))))
  colnames(tbl) <- colnames(obs)
  
  ## Use paste to convert each row into a character string
  
  tblRow <- do.call(paste,c(data.frame(tbl),sep="."))
  obsRow <- do.call(paste,c(data.frame(obs),sep="."))
  
  d <- integer(nrow(tbl)) ## initialize vector of counts
  counts <- table(obsRow) ## Let (the fast) table() do the work
  d[match(names(counts),tblRow)] <- counts ## vector of counts
  out <- list(tbl=as.data.frame(tbl),d=d) ## 'tbl' has to be a dataframe to work with glm()
  class(out) <- "gunter"
  return(out)
}
  
".gunter_MB" <- function(obs){  ## another modification.  Thanks again, Bert!
  # here, 'obs' is an object of class MB.
  f <- function(n){0:n}  # has to start at 0, as 0 count is OK
  tbl <- as.matrix(do.call("expand.grid",sapply(getM(obs),f)))
  cbs <- counts(obs)
  
  ## Use paste to convert each row into a character string
  
  tblRow <- do.call(paste,c(data.frame(tbl),sep="."))
  obsRow <- do.call(paste,c(data.frame(cbs),sep="."))
  
  d <- integer(nrow(tbl)) ## initialize vector of counts
  counts <- table(obsRow) ## Let (the fast) table() do the work
  d[match(names(counts),tblRow)] <- counts ## vector of counts
  out <- list(tbl=as.data.frame(tbl),d=d,m=getM(obs))
  class(out) <- "gunter_MB"
  return(out)
}

".gunter_MB_to_Oarray" <- function(obs){ # convert from an Oarray to a gunter.
  tbl <- which(obs<Inf,arr.ind=TRUE)-1L
  colnames(tbl) <- names(dimnames(obs))

  d <- c(obs)

  m <- dim(obs)-1L
  names(m) <- names(dimnames(obs))

  out <- list(tbl=tbl, d=d, m=m)
  class(out) <- "gunter_MB"
  return(out)
}

".Oarray_to_MB" <- function(a){  #coerces an Oarray object to an MB object
  stopifnot(is.Oarray(a))
  dep <- which(a < Inf,arr.ind=TRUE)-1L
  n <- c(a)
  m <- dim(a)
  names(m) <- names(dimnames(a))
  dep <- dep[rep(seq_along(n),times=n),]
  rownames(dep) <- NULL
  return(MB(dep=dep, m=m, pnames=names(m)))
}
  
setMethod("gunter","data.frame", .gunter)
setMethod("gunter","matrix", .gunter)
setMethod("gunter","MB", .gunter_MB)
setMethod("gunter","Oarray", .gunter_MB_to_Oarray)

"print.gunter" <- function(x, ...){
  print(unclass(x))
  return(invisible(x))
}

"print.gunter_MB" <- function(x, ...){
  print(unclass(x))
  return(invisible(x))
}

"as.array.MB" <- function(x, ...){
  as.array(gunter(x))
}

"as.array.gunter_MB" <- function(x, ...) {  ## "x" is a gunter_MB object, this function coerces x to an Oarray
  m <- x$m
  out <- Oarray(0,m+1,offset=0)
  out[as.matrix(x$tbl)] <- x$d
  jj <- sapply(m,function(i){0:i},simplify=FALSE)
  names(jj) <- names(m)
  dimnames(out) <- jj
  return(out)
}

"Lindsey" <- function(obs,n=NULL,give_fit=FALSE){
  if(inherits(obs,"gunter")){
    return(Recall(obs=obs$tbl,n=obs$d,give_fit=give_fit))
  }
  
  if(!is.null(n)){
    obs <- obs[rep(seq_len(nrow(obs)),n),]
    }
  jj <- gunter(obs)
  Off <- -rowSums(lfactorial(jj$tbl))
  fit <- glm(jj$d ~ -1+offset(Off)+(.)^2,data=data.frame(jj$tbl),family=poisson)
  w <- ncol(obs)
  jj <- exp(fit$coefficients)
  l <- seq_len(w)
  jj[l] <-   jj[l] /sum(jj[l])
  out <- paras(jj[-w])
  pnames(out) <- colnames(obs)
  if(give_fit){
    out <- list(MLE = out, fit = fit)
    class(out) <- "Lindsey_output"
  }
  return(out)
}

"print.Lindsey_output" <- function(x, ...){
  jj <- unclass(x)
  jj[[2]] <- summary(jj[[2]])
  print(jj)
  return(invisible(x))
}
  
"Lindsey_MB" <- function(a){

  if(!inherits(a,"gunter_MB")){
    a <- gunter(a)
  }
  m <- a$m
  if(length(m) != 2){
    stop("only bivariate case coded up")
  }
  o <- names(m)
  d <- a$d
  jj <- a$tbl
  rownames(jj) <- NULL
  x <- data.frame(cbind(d,
                        jj[,1],   # bacon
                        jj[,2],   # eggs
                        jj[,1]*(m[1]-jj[,1]),  # bacon:nbacon
                        jj[,2]*(m[2]-jj[,2]),  # eggs:neggs
                        jj[,1]*jj[,2]          # bacon:eggs
                        ))
  rownames(x) <- NULL
  colnames(x) <-
    c("d",o[1], o[2], paste(o[1],":n",o[1],sep=""), paste(o[2],":n",o[2],sep=""), paste(o[1],":",o[2],sep=""))
  
  Off <- lchoose(m[1],jj[,1]) + lchoose(m[2],jj[,2])
  fit <- glm(d~(.), data=x, family=poisson,offset=Off)
  return(fit)
}

"rMM" <- function(n, Y, paras, burnin=4*Y, every=4*Y, start=NULL){  # Y=100,x
  w <- length(paras)

  update <- function(xin){ # returns an updated RV using xin as a
                           # start; a single iteration of MH
    jj <- sample(w,2,replace=FALSE)
    kernel <- rep(0L,w)
    kernel[jj] <- c(-1L,+1L)
    proposed <- xin + kernel
    if(any(proposed<0)){
      lognum <- -Inf
    } else {
      lognum <- MM_single(proposed, paras,givelog=TRUE)
    }
    logden <- MM_single(xin,paras,givelog=TRUE)
    if ((lognum == -Inf) & (logden == -Inf)) {
      print("this cannot happen")
      alpha <- 0
    } else {
      alpha <- min(1, exp(lognum-logden))
    }
    if (runif(1) < alpha) {
      ans <- proposed
    } else {
      ans <- xin
    }
    return(ans) # NB: this is update() returning
  }

  if(is.null(start)){
    start <- tabulate(sample(w,size=Y,prob=p(paras),replace=TRUE),nbins=w)
  }

  for(i in seq_len(burnin)){
    start <- update(start)
  }

  out <- matrix(0,n,w)
  colnames(out) <- pnames(paras)
  for(i in seq_len(n)){
    for(j in seq_len(every)){
      start <- update(start)
    }
    out[i,] <- start
  }

  return(out)
}
  
"suffstats" <- function(y,n=NULL){
  if(inherits(y,"gunter")){
    return(Recall(y=y$tbl, n=y$d))
  }
  stopifnot(minmax(rowSums(y)))
  if(is.null(n)){n <- rep(1,nrow(y))}
  nc <- ncol(y)
  m <- matrix(0,nc,nc)
  for(i in seq_len(nc)){
    for(j in seq_len(nc)){
      m[i,j] <- sum(y[,i]*y[,j]*n)
    }
  }
  rownames(m) <- colnames(y)
  colnames(m) <- colnames(y)
  out <- 
    list(Y           = sum(y[1,]),
         nobs        = sum(n),
         row_sums    = colSums(sweep(y,1,n,"*")),
         cross_prods = m
         )
  class(out) <- "suffstats"
  return(out)
}

"print.suffstats" <- function(x, ...){
  print(unclass(x))
  return(invisible(x))
}

"summary.suffstats" <- function(object, ...){
  object$row_sums <- object$row_sums/object$nobs
  object$cross_prods <- object$cross_prods/(object$nobs)

  object$Y <- NULL
  object$nobs <- NULL
  class(object) <- "summary.suffstats"

  return(object)
}
 
"expected_suffstats" <- function(L,Y){
  stopifnot(inherits(L,"paras"))
  a <- pnames(L)
  cc <- compositions(Y,length(a))
  rownames(cc) <- a
  probs <- apply(cc,2,MM_single,paras=L)
  probs <- probs/sum(probs)
  rs <- rowSums(sweep(cc,2,probs,"*"))

  f <- function(i){  outer(cc[,i],cc[,i])  }
  
  x <- do.call("abind",c(sapply(seq_len(dim(cc)[2]),f,simplify=FALSE),along=3))
  cp <- apply(sweep(x,3,probs, "*"),1:2,sum)

  return(list(row_sums=rs, cross_prods=cp))
}

"print.summary.suffstats" <- function(x, ...){
  print(unclass(x))
  return(invisible(x))
}
