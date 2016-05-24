
update.rcox <- function(object,
                        vcc       = NULL,
                        ecc       = NULL,
                        splitecc  = NULL,
                        splitvcc  = NULL,
                        joinvcc   = NULL,
                        joinecc   = NULL,
                        addecc    = NULL,
                        dropecc   = NULL,
                        Kstart    = NULL,
                        fit       = TRUE,
                        control   = NULL,
                        trace=object$trace,
                        ...){


  .intersectListList <- function(addccV, oldccV){
    x <- lapply(addccV,
                function(ccnew){
                  lapply(oldccV,
                         function(ccold){ intersect(ccnew,ccold) })})
    x
  }

  nodes <- object$nodes
  
  if (!is.null(joinvcc)){
    oldcc   <- getSlot(object,"vcc")
    joinvcc   <- .ccl2names(joinvcc, oldcc)
    new.ccl   <- .joincc(joinvcc, oldcc)
    new.ccl   <- .addccnames(new.ccl, type="vcc")
    vcc       <- new.ccl
    if (trace>=1)cat(".joining vcc:", toLisp(joinvcc),"\n")
  }

  if (!is.null(joinecc)){
    oldcc   <- getSlot(object,"ecc")
    joinecc   <- .ccl2names(joinecc, oldcc)
    new.ccl   <- .joincc(joinecc, oldcc)
    new.ccl   <- .addccnames(new.ccl, type="ecc")
    ecc       <- new.ccl
    if (trace>=1)cat(".joining ecc:", toLisp(joinecc),"\n")
  }
  
  if (!is.null(splitvcc)){
    oldcc    <- getSlot(object,"vcc")
    splitvcc   <- .ccl2names(splitvcc, oldcc)
    new.ccl    <- .splitcc(splitvcc, oldcc)
    new.ccl    <- .addccnames(new.ccl, type="vcc")
    vcc        <- new.ccl
    if (trace>=1)cat(".splitting vcc:", toLisp(splitvcc),"\n")
  }
  
  if (!is.null(splitecc)){
    oldcc   <- object$ecc
    oldccV  <- object$intRep$eccV
    splitecc   <- .ccl2names(splitecc, oldcc)    
    
    new.ccl    <- .splitcc(splitecc, oldcc)
    new.ccl    <- .addccnames(new.ccl, type="ecc")
    ecc        <- new.ccl
    if (trace>=1)cat(".splitting ecc:", toLisp(splitecc),"\n")
  }

  if (!is.null(addecc)){
    oldcc   <- object$ecc
    oldccV  <- object$intRep$eccV
    addecc  <- formula2names(addecc)
    #cat("addecc:\n"); print(addecc)
    #cat("oldecc:\n"); print(oldcc)
    if (length(oldcc)==0){
      ecc <- addecc
    } else {
      addccV <- indices2vectors(names2indices(formula2names(addecc), nodes, matrix=TRUE))
      x      <- .intersectListList(addccV, oldccV)
      if (length(unlist(x))>0)
        stop("Can not add ecc to model\n");
      ecc <- c(oldcc, addecc)
    } 
    ecc    <- .addccnames(ecc, type="ecc")
    #print(ecc)
    if (trace>=1)cat(".add ecc:", toLisp(addecc),"\n")
  }
  
  if (!is.null(dropecc)){
    oldcc   <- object$ecc
    oldccV  <- object$intRep$eccV
    dropecc  <- formula2names(dropecc)

    dropccV <- indices2vectors(names2indices(formula2names(dropecc), nodes, matrix=TRUE))
    #print(dropccV)
    #print(oldccV)
    
                                        #    idx       <- sapply(dropecc, matchLL2, oldcc)
    idx <- unlistPrim(lapply(dropecc, matchLL2, oldcc))
    idx       <- which(!is.na(idx))
    dropecc   <- dropecc[idx]    
    #idx       <- sapply(dropecc, matchLL2, oldcc)
    idx       <- unlistPrim(lapply(dropecc, matchLL2, oldcc))
    
    new.ccl   <- oldcc[-idx]    
    ecc       <- new.ccl

    if (trace>=1)cat(".drop ecc:", toLisp(dropecc),"\n")
  }

  if (!is.null(vcc)){
    vcc <- .addccnames(vcc, "vcc")
    object$vcc <- vcc
  }

  if (!is.null(ecc)){
    ecc <- .addccnames(ecc, "ecc")
    object$ecc <- ecc
  }
  
  if (trace>=3)
    cat("...(update) Updating internal representation of model object...\n")
  intRep <- .buildInternalRepresentation(vccN      = object$vcc,
                                         eccN      = object$ecc,
                                         dataNames = object$dataRep$dataNames,
                                         trace     = 2)
  object$intRep <- intRep
  object$Kstart <- Kstart
  
  if (!is.null(control)){
    object$control[(namc <- names(control))] <- control
  }
  
  if (fit)# & !is.null(object$fitInfo))
    object$fitInfo <- fit(object, trace=trace, returnModel=FALSE)
  else
    object$fitInfo <- NULL
  return(object)
}








.splitcc <- function(cc, old.ccl){
  idx       <- sapply(cc, matchLL2, old.ccl)
  if ((length(idx)>1) || (!is.na(idx))){
    old.ccl <- old.ccl[-idx]
  }   
  new.cc    <- lapply(unlist(cc, recursive=FALSE),list)
  new.ccl   <- unionL2L2(old.ccl,new.cc)
  new.ccl
}

.joincc <- function(cc, old.ccl){
  idx       <- sapply(cc, matchLL2, old.ccl)
  if ((length(idx)>1) || (!is.na(idx))){
    old.ccl <- old.ccl[-idx]
  }   
  new.cc         <- list(unlist(cc, recursive=FALSE))
  new.ccl        <- unionL2L2(old.ccl, new.cc)
  class(new.ccl) <- union(class(old.ccl),class(cc))
  new.ccl
}

.ccl2names <- function(x,y){
  if (class(x)[1]=="formula")
    x <- list(x)
  cc <- lapply(x, function(xx){
    switch(class(xx),
      "list"    = {xx},
      "numeric" = {y[[xx]]},
      "formula" = {formula2names(xx)}
      )
  })
  cc
}



