#########################################################################
##
## Backward and forward stepwise model selection for iModel's
##
#########################################################################

stepwise.iModel <- function(object,
                          criterion="aic",
                          alpha=NULL,
                          type="decomposable",
                          search="all",
                          steps = 1000,
                          k = 2,
                          direction="backward",
                          fixinMAT=NULL,
                          fixoutMAT=NULL,
                          details=0,
                          trace=2, ...)
{

  direction  <- match.arg(direction, c("backward", "forward", "both"))
  criterion  <- match.arg(criterion, c("aic",        "test"))
  search     <- match.arg(search,    c("headlong",   "all"))

  if (isGSD_glist(object$glist)[2]){
    type <- match.arg(type,      c("decomposable", "unrestricted"))
  } else {
    type <- "unrestricted"
  }

  if(details>=1){
    cat("STEPWISE: ",
        "\n criterion:", criterion, if (criterion=="aic") paste("( k =", round(k,2),")"),
        "\n direction:", direction,
        "\n type     :", type,
        "\n search   :", search,
        "\n steps    :", steps, "\n")
  }

  if (direction=="backward"){
    ans <- backward(object, type=type, search=search,
                    criterion=criterion, alpha=alpha, steps=steps,
                    k=k, fixinMAT=fixinMAT, details=details, trace=trace)#$object
  } else {
    ans <- forward(object, type=type,  search=search,
                   criterion=criterion, alpha=alpha, steps=steps,
                   k=k, fixoutMAT=fixoutMAT, details=details, trace=trace)#$object
  }
  return(ans)
}



backward <- function(object, criterion="aic", alpha=NULL, type="decomposable", search="all",
                     steps=1000,  k=2, fixinMAT=NULL, details=1, trace=2, ...)
{
  type   <- match.arg(type,   c("decomposable", "unrestricted"))
  search <- match.arg(search, c("headlong",     "all"))

  ## Need this because checking for decomposability is special for mixed models.
  if (inherits(object, "mModel")){
    discrete <- object$datainfo$disc.names  
  } else {
    discrete <- NULL
  }
  
  isgsd  <- isGSD_glist(object$glist, discrete=discrete)
  vn     <- object$varNames

  if (!all(isgsd)){
    type <- "unrestricted"
  }

  if (is.null(alpha)){
    alpha <- if (criterion=="aic") 0 else 0.05
  }

  switch(criterion,
       "aic" ={opt.op    <- which.min
               comp.op   <- `<`
               outstring <- "change.AIC"
               crit.str  <- "aic"},
       "test"={opt.op    <- which.max
               comp.op   <- `>`
               outstring <- "p.value"
               crit.str  <- "p.value"
               })

  itcount <- 1
  t0 <- proc.time()

  testFun <- switch(search,
                     "headlong" = {.testInEdges_headlong},
                     "all"      = {.testInEdges_all})
  
  .infoPrint(details, 1,
             if (criterion=="aic"){
               cat(sprintf("BACKWARD: type=%s search=%s, criterion=%s(%4.2f), alpha=%4.2f \n",
                           type, search, criterion, k, alpha ))
             } else {
               cat(sprintf("BACKWARD: type=%s search=%s, criterion=%s, alpha=%4.2f \n",
                           type, search, criterion, alpha ))
             })
  .infoPrint(details, 1,
             cat(sprintf("Initial model: is graphical=%s is decomposable=%s\n",
                         isgsd[1], isgsd[2])))

  amat    <- glist2adjMAT(object$glist)
  edgeMAT <- getEdges(amat, type=type, ingraph=TRUE, discrete=discrete)
  fixNUM  <- pairs2num(fixinMAT,  vn)
  edgeMAT <- .subtract.fix(fixNUM,  edgeMAT, vn)

  if (nrow(edgeMAT)==0){
    if (details>=1)
      cat(sprintf("No edges can be removed\n"))
    return(NULL)
  }
  
  repeat{
    testMAT <- testFun(object, edgeMAT, comp.op=comp.op, crit.str=crit.str,
                       alpha=alpha, k=k, amat=amat, ...)
    if (details>=2) print(testMAT, row.names=FALSE, digits=4)
    
    statvec   <- testMAT[,crit.str]
    opt.idx   <- opt.op(statvec)
    
    if (comp.op( statvec[opt.idx], alpha)) {
      opt.edge    <- as.character(testMAT[opt.idx,c("V1","V2")])
      
      if (details>=1) cat(sprintf("  %s %9.4f Edge deleted: %s\n",
            outstring, statvec[opt.idx], .toString(opt.edge)))

      ## Update model object and prepare to start all over again:
      object  <- update(object, list(drop.edge=opt.edge))

      amat    <- glist2adjMAT(object$glist)
      edgeMAT <- getEdges(amat, type=type, ingraph=TRUE, discrete=discrete)
      edgeMAT <- .subtract.fix(fixNUM, edgeMAT, vn)
      if (nrow(edgeMAT)==0 | itcount==steps) { break }
    } else {
      break
    }
    itcount <- itcount + 1
  }
##  cat(sprintf("itcount %i\n", itcount))
  return(object)
}



forward <- function(object, criterion="aic", alpha=NULL, type="decomposable", search="all",
                     steps=1000,  k=2, fixoutMAT=NULL, details=1, trace=2, ...)
{
  type   <- match.arg(type,   c("decomposable", "unrestricted"))
  search <- match.arg(search, c("headlong",     "all"))

  ## Need this because checking for decomposability is special for mixed models.
  if (inherits(object, "mModel")){
    discrete <- object$datainfo$disc.names  
  } else {
    discrete <- NULL
  }
  
  isgsd  <- isGSD_glist(object$glist, discrete=discrete)
  vn     <- object$varNames

  if (!all(isgsd)){
    type <- "unrestricted"
  }

  if (is.null(alpha)){
    alpha <- if (criterion=="aic") 0 else 0.05
  }

  switch(criterion,
         "aic" ={opt.op    <- which.min
                 comp.op   <- `<`
                 outstring <- "change.AIC"
                 crit.str  <- "aic"},
         "test"={opt.op    <- which.min
                 comp.op   <- `<`
                 outstring <- "p.value"
                 crit.str  <- "p.value"
               })

  itcount <- 1
  t0 <- proc.time()

  testFun <- switch(search,
                     "headlong" ={.testOutEdges_headlong},
                     "all"      ={.testOutEdges_all})

  .infoPrint(details, 1,
             if (criterion=="aic"){
               cat(sprintf("FORWARD: type=%s search=%s, criterion=%s(%4.2f), alpha=%4.2f \n",
                           type, search, criterion, k, alpha ))
             } else {
               cat(sprintf("FORWARD: type=%s search=%s, criterion=%s, alpha=%4.2f \n",
                           type, search, criterion, alpha ))
             }
             )
  .infoPrint(details, 1,
             cat(sprintf("Initial model: is graphical=%s is decomposable=%s\n",
                         isgsd[1], isgsd[2])))

  amat    <- glist2adjMAT(object$glist)
  edgeMAT <- getEdges(amat, type=type, ingraph=FALSE, discrete=discrete)
  fixNUM  <- pairs2num(fixoutMAT,  vn)
  edgeMAT <- .subtract.fix(fixNUM,  edgeMAT, vn)

  if (nrow(edgeMAT)==0){
    if (details>=1)
      cat(sprintf("No edges can be added\n"))
    return(NULL)
  }
  
  repeat{
    testMAT   <- testFun(object, edgeMAT, comp.op=comp.op, crit.str=crit.str,
                         alpha=alpha, k=k, amat=amat, ...)
    if (details>=2) print(testMAT,row.names=FALSE, digits=4)
    
    statvec   <- testMAT[,crit.str]
    opt.idx   <- opt.op(statvec)

    if (comp.op( statvec[opt.idx], alpha)) {
      opt.edge    <- as.character(testMAT[opt.idx,c("V1","V2")])
      
      if (details>=1) cat(sprintf("  %s %9.4f Edge added: %s\n",
            outstring, statvec[opt.idx], .toString(opt.edge)))

      ## Update model object and prepare to start all over again:
      object  <- update(object, list(add.edge=opt.edge))
      amat    <- glist2adjMAT(object$glist)
      edgeMAT <- getEdges(amat, type=type, ingraph=FALSE, discrete=discrete)
      edgeMAT <- .subtract.fix(fixNUM, edgeMAT, vn)
      if (nrow(edgeMAT)==0 | itcount==steps) { break }
    } else {
      break
    }
    itcount <- itcount + 1
  }
  return(object)
}



## Used for modifying an edgeMAT by subtracting those edges in fixMAT from edgeMAT
.subtract.fix <- function(fixNUM, edgeMAT, vn){

    if (is.null(fixNUM)){
        return(edgeMAT)
    } else {
        if (nrow(edgeMAT)==0){
            return(edgeMAT)
        } else {
                                        #fixNUM  <- pairs2num(fixMAT,  vn)
                                        #print(fixNUM)
          edgeNUM <- pairs2num(edgeMAT, vn)
                                        #print(edgeNUM)
          edgeMAT[(1:nrow(edgeMAT))[-na.omit(matchPrim(fixNUM,edgeNUM))],,drop=FALSE]
        }
    }
}




