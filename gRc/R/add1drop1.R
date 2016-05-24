## Calculate test statistic for addition of each atomic ecc not already in the model
## Note: Based on fitting a new model for each new ecc, and may hence be slow...
##

add1.rcox <- function(object, scope, details=0, trace=0, ...){
  ##cat("add1.rcox\n")
  n       <- dataRep(object,"n")
  if (missing(scope)){
    eNew <- getedges(object, complement=TRUE)
  } else{
    eNew <- .addccnames(formula2names(scope),type="ecc")
  }
  
  if (is.L(eNew)){
    eNew <- lapply(eNew, list)
    class(eNew)<- c("colourClass", "list")
  }
  
  if (length(eNew)==0)
    return(NULL)

  object$control$vcov <- NULL

  res <- rep(NA, length(eNew))
  for (i in 1:length(eNew)){
    e      <- eNew[i]
    mtmp   <- update(object, addecc=e)
    dev    <- -2*(logL(object)-logL(mtmp))
    res[i] <- dev
  }

  eNew <- .addccnames(eNew,"ecc")
  ans  <- data.frame(cc=names(eNew), X2=res, df=1)
  ans  <- .addStat(ans, n=n, direction="add")
  
  attr(ans,"ccterms") <- eNew
  ans2 <- structure(list(tab=ans, cc=eNew, details=details),
                    class=c("statTable","data.frame"))
  ans2
}


.compareModelRCOX <- function(x, basestat=c(0,0), criterion="aic", k=2, direction="drop"){

  criterion <- match.arg(criterion, c("aic","test")) 
  direction <- match.arg(direction, c("drop","add"))

  switch(criterion,
    "aic"={
      if (direction=="drop")
        ans <- c(-x[1]+basestat[1] + k*(x[2]-basestat[2]))
      else
        ans <- c(-x[1]+basestat[1] + k*(x[2]-basestat[2]))
  

    },
    "test"={
      ans <- 1-pchisq(abs(x[1]-basestat[1]), df=abs(x[2]-basestat[2]))
    } 
  )
  return(ans)
}























print.statTable <- function(x,...){
  print(x$tab)
  if (x$details>=1){
    if (!is.null(x$cc)){
      cat("\ncc:\n")
      print(x$cc)
    }
    if (!is.null(x$cc1)){
      cat("\ncc1:\n")
      print(x$cc1)
      cat("cc2:\n")
      print(x$cc2)
    }
  }
  cat(paste("\nAvailable components:", paste(setdiff(names(x),"details"),collapse=' ')),"\n")
  return(invisible(x))
}
  

## Calculate test statistic for deletion of each ecc in the model
## Note: By default based on Wald statistics from existing model
##
drop1.rcox <- function(object, scope, details=0, trace=0, stat="wald", ...){
  stat <- match.arg(stat,c("wald","dev"))
  n    <- dataRep(object,"n")

  if (missing(scope))
    ec  <- getSlot(object,'ecc')
  else
    ec  <- .addccnames(formula2names(scope),type="ecc")
  
  if (details>=1){
    cat("Statistic:", stat, "\n")
  }
  
  if (length(ec)==0)
    return(NULL)

  res <- rep(NA, length(ec))

  if (stat=="wald"){
    V   <- vcov(object)
    b   <- coef(object)
    ofs <- length(getSlot(object,"vcc"))
  
    ccidx <- sapply(ec, matchLL2, ec) # FIXME: What is this???
    lcc   <- length(ccidx)
    for (i in ccidx){    
      i2   <- i + ofs
      dev  <- b[i2]^2/V[i2,i2]
      res[i] <- dev
    }
  } else {
    for (i in 1:length(ec)){
      e      <- ec[[i]]
      mtmp   <- update(object, dropecc=list(e))
      dev    <- 2*(logL(object)-logL(mtmp))
      res[i] <- dev
    }
  }
  
  ans <- data.frame(cc=names(ec), X2=res, df=1)
  ans <- .addStat(ans, n=n, direction="drop")
  attr(ans,"ccterms") <- ec
  ans
  ec <- .addccnames(ec,"ecc")
  ans2 <- structure(list(tab=ans, cc=ec, details=details), class=c("statTable","data.frame"))
  ans2
}


## Evaluate edge *not in* the model
##
evalOutECC <- function(object, edlist,
                         alpha=0, criterion="aic", k=2,
                         headlong=FALSE,
                         random=TRUE,
                         print=1, ...){

  ##print("evalOutEdges")
  ##bstat     <- extractFIT(object$fitinfo)
  statlist  <- list()
  optModel  <- object
  optStat   <- 999999
  optEdge   <- NULL

  criterion <- match.arg(criterion, c("aic","test"))
  if (headlong && random){
    edlist <- edlist[sample(length(edlist), replace=FALSE)]
  }

  ##cat("bstat:\n"); print(bstat)
  res <- rep(NA, length(edlist))
  for (ii in seq_along(edlist)){
    ed     <- edlist[[ii]]
    e      <- edlist[ii]
    mtmp   <- update(object, addecc=e)
    dev    <- -2*(logL(object)-logL(mtmp))
    res[ii] <- dev
    
    stat   <- .compareModelRCOX(c(dev,1), c(0,0),criterion=criterion, k=k,direction="add")

    if (print>=1){
      sss <- sprintf(" %30s %15.10f %s", paste(unlist(e),collapse=" ~ "),
                     stat, ifelse(stat<alpha,"+","-"))
      cat(sss, "\n")
    }
    if (stat<min(optStat,alpha)){ ## min because we minimize
      optModel <- mtmp
      optStat  <- stat
      optEdge  <- ed
    }

    if (headlong && optStat<alpha) ## < because we minimize
      break()
    statlist[[ii]] <- stat
  }

  ans <- list(optEdge=optEdge, optStat=optStat, optModel=optModel,
              statlist=unlist(statlist))
  
  return(ans)
}


## Evaluate edge *in* the model
##
evalInECC <- function(object, edlist,
                        alpha=0, criterion="aic", k=2,
                        stat="wald",
                        headlong=FALSE,
                        random=TRUE,
                        print=1, ...){

  criterion <- match.arg(criterion, c("aic","test"))
  stat <- match.arg(stat,c("wald","dev"))

  if (headlong && random){
    edlist <- edlist[sample(length(edlist), replace=FALSE)]
  }

  
  statlist  <- list()
  optModel  <- object
  optStat   <- -9999
  optEdge   <- NULL

  res <- rep(NA, length(edlist))

  if (stat=="wald"){
    V   <- vcov(object)
    b   <- coef(object)
    ofs <- length(getSlot(object,"vcc"))
    ccidx <- sapply(edlist, matchLL2, object$ecc) # Needed for wald with random
    lcc   <- length(ccidx)
  }


  for (ii in 1:length(edlist)){
    ed      <- edlist[[ii]]
    if (stat=="wald"){
      i    <- ccidx[ii]
      i2   <- i + ofs
      dev  <- b[i2]^2/V[i2,i2]
      res[i] <- dev      
    } else{
      mtmp   <- update(object, dropecc=list(ed))
      dev    <- 2*(logL(object)-logL(mtmp))
      res[ii] <- dev
    }


    statValue <- .compareModelRCOX(c(dev,1), c(0,0),criterion=criterion, k=k,...)
    
    if (print>=1){
      sss <- sprintf(" %30s %10.4f %15.10f %s", paste(unlist(ed),collapse=" ~ "),
                     dev, statValue, ifelse(statValue>alpha,"+","-"))
      cat(sss, "\n")
    }


    if (statValue>max(optStat,alpha)){ ## max because we maximize
      optStat  <- statValue
      optEdge  <- ed
      if (stat=="dev")
        optModel <- mtmp
    }

    statlist[[ii]] <- statValue
    if (headlong && optStat>alpha) ## > because we maximize
      break()

  }

  if (stat=="wald"){
    if (length(optEdge))
      optModel <- update(object, dropecc=list(optEdge))
    else
      optModel <- object
  }   

  ans <- list(optEdge=optEdge, optStat=optStat, optModel=optModel,
              statlist=unlist(statlist))
  return(ans)
}




## ### PORTO - working on this one
## ###
## add1.rcox <- function(object, scope, details=0, trace=0, ...){
##   ##cat("add1.rcox\n")
##   n       <- dataRep(object,"n")
##   if (missing(scope)){
##     eNew <- getedges(object, complement=TRUE)
##   } else{
##     eNew <- .addccnames(formula2names(scope),type="ecc")
##   }
  
##   if (is.L(eNew)){
##     eNew <- lapply(eNew, list)
##     class(eNew)<- c("colourClass", "list")
##   }
  
##   if (length(eNew)==0)
##     return(NULL)

##   object$control$vcov <- NULL

##   evalOutEdges(object, eNew, crit="test",alpha=0.05)

##   ## res <- rep(NA, length(eNew))
## ##   for (i in 1:length(eNew)){
## ##     e      <- eNew[i]
## ##     mtmp   <- update(object, addecc=e)
## ##     dev    <- -2*(logL(object)-logL(mtmp))
## ##     res[i] <- dev
## ##   }

## ##   eNew <- .addccnames(eNew,"ecc")
## ##   ans  <- data.frame(cc=names(eNew), X2=res, df=1)
## ##   ans  <- .addStat(ans, n=n, direction="add")
  
## ##   attr(ans,"ccterms") <- eNew
## ##   ans2 <- structure(list(tab=ans, cc=eNew, details=details),
## ##                     class=c("statTable","data.frame"))
## ##   ans2
## }

