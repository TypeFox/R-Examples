##### Distributed Logistic Multinomial Regression  ######

## define class
setClass("dmrcoef", contains="dgCMatrix")

## inner loop function
onerun <- function(xj, argl){
  argl$y <- xj
  if(argl$nzcheck | argl$mlcheck){ 
    xnz <- as.matrix(argl$x[drop(argl$y>0),,drop=FALSE])
    if(nrow(xnz)==0) return(NULL)
    if(argl$nzcheck) fullrank <- which(colSums(xnz!=0)!=0)
    else{
      Q <- qr(cbind(1,xnz))
      fullrank <- Q$pivot[ 2:Q$rank ] - 1
    }
    if(is.null(argl$varweight)) argl$varweight <- rep(1,ncol(argl$x))
    argl$varweight[-fullrank] <- Inf
    argl$free <- argl$free[ (argl$free %in% fullrank) ]
  }
  if(argl$cv) fit <- do.call(cv.gamlr,argl)
  else fit <- do.call(gamlr,argl)
  fit$nobs <- argl$nobs
  return(fit)
}

## main function
dmr <- function(cl, covars, counts, mu=NULL, bins=NULL, verb=0, cv=FALSE, ...)
{
  if(!is.null(cl)){
    if(!inherits(cl,"cluster")) stop("first argument `cl' must be NULL or a socket cluster.")
  }
  #build the default argument list
  argl <- list(...)
  argl$family <- "poisson"
  if(is.null(argl$nlambda))
    argl$nlambda <- formals(gamlr)$nlambda
  argl$verb <- max(verb-1,0)
  argl$cv <- cv
  if(is.null(argl$mlcheck))
    argl$mlcheck <- FALSE
  if(is.null(argl$nzcheck))
    argl$nzcheck <- FALSE
  ## collapse and clean
  chk <- collapse(covars, counts, mu, bins)
  if(verb)
    cat(sprintf("fitting %d observations on %d categories, %d covariates.\n",
        nrow(chk$v), ncol(chk$counts), ncol(chk$v)))
  argl$x <- chk$v
  argl$shift <- chk$mu
  argl$nobs <- sum(chk$nbin)
  p <- ncol(chk$counts)
  vars <- colnames(chk$counts)
  ## cleanup
  rownames(argl$x) <- rownames(chk$counts) <- NULL
  counts <- chk$counts
  rm(covars,mu,chk)

  ## convert X to list
  if(verb) cat("converting counts matrix to column list...\n")
  C <- ifelse(is.null(cl),Inf,length(cl))
  if(C < p/4){
    chunks <- round(seq(0,p,length.out=C+1))
    counts <- lapply(1:C, 
      function(i) counts[,(chunks[i]+1):chunks[i+1]])
    counts <- parLapply(cl,
                counts, 
                function(x) 
                  sapply(colnames(x), 
                  function(j) x[,j,drop=FALSE]))
    counts <- unlist(counts,recursive=FALSE)
  } else{
    counts <- sapply(vars,
      function(j) counts[,j,drop=FALSE]) }

  ## lapply somehow, depending on cl and p
  if(is.null(cl)){
    if(verb) cat("running in serial.\n ")
    mods <- lapply(counts,onerun,argl=argl) 
  }
  else{
    if(verb){ 
     cat("distributed run.\n") 
     print(cl) }
    mods <- parLapply(cl,counts,onerun,argl=argl) 
  }
    
  ## classy exit
  class(mods) <- "dmr"
  attr(mods,"nobs") <- argl$nobs
  attr(mods,"nlambda") <- argl$nlambda
  attr(mods,"mu") <- argl$shift
  return(mods)
}

coef.dmr <- function(object, ...){
  B <- lapply(object,coef, ...)
  B[sapply(B,is.null)] <- Matrix(0)
  bx <- unlist(lapply(B,function(b) b@x))
  bi <- unlist(lapply(B,function(b) b@i))
  bp <- c(0,
    cumsum(unlist(lapply(B,function(b) b@p[-1]))))
  Bs <- sparseMatrix(i=bi+1,p=bp,x=bx,
    dims=c(nrow(B[[1]]),length(B)),
    dimnames=list(rownames(B[[1]]),names(B)))
  Bs <- as(as(Bs,"dgCMatrix"),"dmrcoef")
  return(Bs)
}

## method predict functions
predict.dmr <- function(object, newdata, 
                  type=c("link","response","class"), ...){
  B <- coef(object, ...)
  predict(B,newdata=newdata,type=type)
}

predict.dmrcoef <- function(object, newdata, 
                  type=c("link","response","class"), ...)
{
  if(inherits(newdata,"simple_triplet_matrix"))
    newdata <- sparseMatrix(i=newdata$i,j=newdata$j,x=newdata$v,
      dims=dim(newdata),dimnames=dimnames(newdata))
  if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
  if(is.data.frame(newdata)){ newdata <- as.matrix(newdata) }

  type=match.arg(type)
  if(type=="reduction")
    stop("type `reduction' has been replaced by the `srproj' function in the textir library.")

  eta <- t(tcrossprod(t(object[-1,,drop=FALSE]),newdata) + object[1,])
  if(type=="response"){
    expeta <- exp(eta)
    eta <- expeta/rowSums(expeta) }
  rownames(eta) <- rownames(newdata)
  colnames(eta) <- colnames(object)
  
  if(type=="class"){
    c <- apply(eta,1,function(e) colnames(eta)[which.max(e)])
    return(c)
  }
  else return(as.matrix(eta))

}

setGeneric("predict")
setMethod("predict","dmrcoef",predict.dmrcoef)
