
## Stepwise joining of existing colour classes in model
##
stepjoin1 <- function(object, scope, type='ecc', criterion='aic', steps=1000, k=2,
                      alpha=0.05, stat="wald",
                      details=1,trace=0)
{
  criterion <- match.arg(criterion, c("test","aic","bic"))
  type      <- match.arg(type, c("ecc","vcc"))
  stat      <- match.arg(stat, c("wald","dev"))

  if (criterion %in% c("aic","bic")){
    alpha <- 0
  }
  
  if (details>=1){   
    typeStr <- ifelse (type=="vcc", "vertex colour classes", "edge colour classes")
    cat("Stepwise joining of", typeStr, "\n")
    cat(.stepDetailsString(criterion, stat, k, alpha),"\n")
  }   

  if (missing(scope)){
    scopecc <- getcc(object,type)
  } else{
    scopecc <- .addccnames(formula2names(scope),type=type)
  }

  currcc <- getcc(object,type)
##  print("Start - scopecc:");  print(scopecc); print("--------------")
  
  x<-sapply(scopecc, function(e){
    i <- match.containsLL2(e, currcc);
    i
  })
  names(x)<-NULL
  ##print(x)
  
  untouch     <- currcc[-x]
  len.untouch <- length(untouch)
  ##print(.addccnames(untouch,type)); print("------------")
  
  changelist  <- NULL
  modelChange <- FALSE



  stepcount   <- 1
  repeat {
    #print("Iteration - scopecc (before)"); print(scopecc)
    cctab <- comparecc(object, cc1=scopecc, cc2=scopecc, type=type, stat=stat,
                       details=details-1)
    #print(cctab)
    #print(cctab$cc1)
    ##print(cctab$cc2)
    tab   <- cctab$tab

    if (is.null(tab))
      return(object)
    
    statvalue <- .extractStat(tab, criterion, k)
    idx       <- which.max(statvalue)
    optstat   <- statvalue[idx]
    if (details>=3)
      print(tab)

    if (optstat > alpha){
      modelChange <- TRUE
      ##cc          <- attributes(tab)$ccterms[[idx]]

      ccid <- as.character(tab[idx,c("cc1","cc2")])
      #print(ccid)
      cc   <- cctab$cc1[ccid]
      
      if (details>=1){
        cat("Joining", paste(ccid, collapse=' '),":",
            paste(formula2string(names2formula(cc)), collapse="; "),"\n")
        xxx <- tab[idx,-c(1:2)]
        cat(.printTestResult(xxx),"\n")
      }
      ##print(cc); print(cctab$cc1[names(cc)]); ##print(cctab$cc2[idx]); 
      
      KS <- fitInfo(object,"K")

      scopecc[names(cc)] <- NULL;      
      scopecc <- .addccnames(.joincc(list(unlist(cc,recursive=FALSE)), scopecc), type)
      #print("Iteration - scopecc (after)"); print(scopecc)

      newcc <- .addccnames(unionL2L2(untouch, scopecc),type)
      #print("newcc"); print(newcc)
      
      switch(type,
             "ecc"={
               object  <- update(object, ecc=newcc, trace=trace, fit=FALSE)
             },
             "vcc"={
               object  <- update(object, vcc=newcc, trace=trace, fit=FALSE)
             })




      
      KS     <- findKinModel(object, KS, type=object$type,regularize=TRUE)      
      object <- update(object, Kstart=KS, fit=TRUE)
      
      stepcount  <- stepcount + 1
      changelist <- c(changelist, cc)

      #currcc <- getcc(object,type)
      ##print("Iterate:");  print(currcc)
      
      if (stepcount>steps)
        break()
    } else {
      if (details>=3){
        print(tab)
      }
      break()
    }
    
  }
  object$change <- changelist
  if (details>=1)
    cat("\n")
  if (modelChange)
    return(object)
}






##
## Stepwise addition of atomic ecc's
##

## stepadd1 <- function(object, criterion="aic", steps=1000, k=2, alpha=0.05, details=1, trace=0){

##   criterion <- match.arg(criterion, c("test","aic","bic")) 
##   if (criterion %in% c("aic","bic")){
##     alpha <- 0
##   }
##   if (details>=1){

##     cat("Stepwise addition of atomic edge colour classes; ")  
##     switch(criterion,
##            "aic"={
##              if (k==2) 
##                cat(paste("criterion=aic",  sep=''), "\n")  
##              else
##                cat(paste("criterion=aic", "k=", k, sep=''), "\n")  
##            },
##            "bic"={
##              cat(paste("criterion=bic", sep=''), "\n")      
##            }, 
##            "test"={
##              cat(paste("criterion=test", " alpha=", alpha,sep=''), "\n")      
##            })

##   }


  
##   modelChange <- FALSE
##   stepcount <- 1
##   repeat{
##     ##tab  <- add1(object)

##     add1obj <- add1(object)

##     print(add1obj)
##     tab     <- add1obj$tab

##     if (is.null(tab))
##       break()
##     stat <- .extractStat(tab, criterion, k)
##     idx  <- which.min(stat)
##     print(idx)
##     optstat <- stat[idx]
##     if (details>=3)
##       print(tab)

##     print(optstat)
##     if (optstat < alpha){
##       print("yes...")
##       ##cc <- attributes(tab)$ccterms[idx]
##       modelChange <- TRUE
##       cc <- add1obj$cc[idx]
##       if (details>=2){
##         #cat("Adding: ", toLisp(cc),"\n\n")
##         cat("Adding:", formula2string(names2formula(cc))[[1]],"\n")
##         ##print(tab[idx,])
##         cat(.printTestResult(tab[idx,-1]),"\n\n")
##       }

##       ##object  <- update(object, addecc=list(cc))
##       object  <- update(object, addecc=cc)
##       stepcount <- stepcount + 1
##       if (stepcount>steps)
##         break()
##     } else {      
##       if (details>=2){
##         print(tab)
##       }
##       break()
##     }
##   }
##   if (details>=1)
##     cat("\n")

  
##   if (modelChange)
##     return(object)
##   #else
##   #  return(NULL)
## }

## Add edge
## Test: add edge with smallest p-value
##    if (stat>max(optStat,alpha)){
##    if (stat>max(optStat,alpha)){




## ## Stepwise deletion of ecc's from the model
## ##
## stepdrop1 <- function(object, criterion='aic', steps=1000, k=2,   alpha=0.05, stat="wald", details=1, trace=0){

##   criterion <- match.arg(criterion, c("test","aic","bic"))
##   stat      <- match.arg(stat, c("wald", "dev"))
##   if (criterion %in% c("aic","bic")){
##     alpha <- 0
##   }
##   if (details>=1){

    
##     #typeStr <- ifelse (type=="vcc", "vertex colour classes", "edge colour classes")
##     cat("Stepwise dropping of edge colour classes\n")
##     switch(criterion,
##            "aic"={
##              if (k==2) 
##                cat(paste("statistic=",stat," criterion=aic",  sep=''), "\n")  
##              else
##                cat(paste("statistic=",stat," criterion=aic", "k=", k, sep=''), "\n")  
##            },
##            "bic"={
##              cat(paste("statistic=",stat," criterion=bic", sep=''), "\n")      
##            }, 
##            "test"={
##              cat(paste("statistic=",stat," criterion=test", " alpha=", alpha, sep=''), "\n")      
##            })
##   }


  
##   changelist <- NULL
##   modelChange <- FALSE
##   stepcount <- 1
##   repeat{

##     drop1obj <- drop1(object, stat=stat, details=details-1)
##     tab <- drop1obj$tab
    
##     if (is.null(tab))
##       break()
##     statvalue <- .extractStat(tab, criterion, k)
##     idx  <- which.max(statvalue)
##     optstat <- statvalue[idx]

##     if (details>=2)
##       print(tab)

##     if (optstat>alpha){
##       modelChange <- TRUE
##       cc          <- drop1obj$cc[idx]
      
##       if (details>=1){
##         cat("Dropping:", formula2string(names2formula(cc))[[1]],"\n")        
##         cat(.printTestResult(tab[idx,-1]),"\n\n")
##       }
##       object     <- update(object, dropecc=cc)
##       changelist <- c(changelist, cc)
##       stepcount  <- stepcount + 1
##       if (stepcount>steps)
##         break()
##     } else {
##       if (details>=2){
##         print(tab)
##       }
##       break()
##     }
##   }
##   object$change <- changelist
##   if (details>=1)
##     cat("\n")

##   if (modelChange)
##     return(object)
## }








## Stepwise splitting of colour classes in model
##
stepsplit1 <- function(object, type='ecc', criterion='aic', steps=1000, k=2, alpha=0.05, stat="wald", details=1, trace=0){

  criterion <- match.arg(criterion, c("test","aic","bic"))
  type      <- match.arg(type, c("ecc","vcc"))
  stat      <- match.arg(stat, c("wald","dev"))

  if (criterion %in% c("aic","bic")){
    alpha <- 0
  }
  
  if (details>=1){
    ###cat("Stepwise splitting of colour classes, type:", type, "\n")  

    typeStr <- ifelse (type=="vcc", "vertex colour classes", "edge colour classes")
    cat("Stepwise splitting of", typeStr, "\n")
    switch(criterion,
           "aic"={
             if (k==2) 
               cat(paste("statistic=",stat," criterion=aic",  sep=''), "\n")  
             else
               cat(paste("statistic=",stat," criterion=aic", "k=", k, sep=''), "\n")  
           },
           "bic"={
             cat(paste("statistic=",stat," criterion=bic", sep=''), "\n")      
           }, 
           "test"={
             cat(paste("statistic=",stat," criterion=test", " alpha=", alpha,sep=''), "\n")      
           })
  }
  
  changelist  <- NULL
  modelChange <- FALSE
  stepcount <- 1
  repeat {
    ccl <- getSlot(object, type)
    splittab <- split1(object, type=type)##  stat=stat)

    tab <- splittab$tab
    
    if (is.null(tab))
      return(object)
    statvalue <- .extractStat(tab, criterion, k)
    idx  <- which.min(statvalue)
    optstat <- statvalue[idx]
    if (details>=3)
      print(tab)


    if (optstat<alpha){
      modelChange <- TRUE
      ##cc <- attributes(tab)$ccterms[[idx]]
      cc <- attributes(tab)$ccterms[idx]
      #cc <<- cc
      #print(cc)
      if (details>=1){
        #cat("Split: ", toLisp(cc),"\n")
        cat("Splitting:",paste(formula2string(names2formula(cc)), collapse="; "),"\n")
      }
      if (details>=1){
        ##print(tab[idx,])
        xxx <- tab[idx,-1]
        cat(.printTestResult(xxx), "\n\n")
      }
      
      switch(type,
             "ecc"={
               ##object  <- update(object, splitecc=list(cc), trace=0)
               object  <- update(object, splitecc=cc, trace=0)
             },
             "vcc"={
               ##object  <- update(object, splitvcc=list(cc), trace=0)
               object  <- update(object, splitvcc=cc, trace=0)
             })


      changelist <- c(changelist, cc)
      stepcount <- stepcount + 1
      if (stepcount>steps)
        break()
    } else {
      if (details>=3){
        print(tab)
      }
      break()
    }
  }
  object$change <- changelist
  if (modelChange)
    return(object)
}





### INTERNALS
###
.printTestResult <- function(xxx){
  xxx<-as.numeric(xxx)
  sprintf("  X2: %f df: %d p: %f aic: %f bic: %f", xxx[1], xxx[2], xxx[3], xxx[4], xxx[5])
}







.stepDetailsString <- function(criterion, stat, k, alpha){
  switch(criterion,
         "aic"={
           if (k==2) 
             (paste("statistic=",stat," criterion=aic",  sep=''))  
           else
             (paste("statistic=",stat," criterion=aic", "k=", k, sep=''))  
         },
         "bic"={
           (paste("statistic=",stat," criterion=bic", sep=''))      
         }, 
         "test"={
           (paste("statistic=",stat," criterion=test", " alpha=", alpha,sep=''))      
         })
}




.extractStat <- function(tab, criterion=c("aic","test"), k=2 , direction="drop"){
  criterion <- match.arg(criterion, c("aic","test")) 
  direction <- match.arg(direction, c("drop","add"))

  #direction=drop: AIC is a maximization problem
  #direction=add : AIC is a minimization problem

  # Note dev=2logL(large)-2logL(small) where large is 'fixed'
  # Maximization: we maximize 2logL(j)+k df(j); that is we must
  # maximize -dev+k df(j)
  switch(criterion,
         "aic"  = {
           if (direction=="drop")
             -tab$X2 + k*tab$df
           else
             -tab$X2 - k*tab$df  ## FIXME : signs???
         },
         "test" = {1-pchisq(tab$X2, tab$df)}
         )  
}




stepadd1 <- function(object,  criterion="aic", steps=1000, k=2, alpha=0.05,
                    headlong=FALSE, random=TRUE, details=1, trace=0,...){

  edlist <- lapply(getedges(object, complement=TRUE), list)
  
  itcount <- 1
  
  t0 <- proc.time()
  repeat{
    ev <- evalOutECC(object, edlist=edlist, headlong=headlong, random=random,
                       alpha=alpha,
                       criterion=criterion, k=k, print=details-1)
    if (is.null(ev$optEdge) || itcount>steps)
      break()

    sss <- sprintf("Added edge: %30s   statistic %10.5f",
                   paste(unlist(ev$optEdge), collapse=" ~ "), ev$optStat)
    cat(sss, "\n")
    object  <- ev$optModel     
    edlist  <- getedges(object, complement=TRUE)
    edlist  <- lapply(edlist, list)
    ###edlist  <- edlist[sample(length(edlist))]
    itcount <- itcount+1
  }
  if (details>=1)
    cat("Time (stepwise forward): ", (proc.time()-t0)[3], "\n")
  
  return(object)
}


## Stepwise deletion of ecc's from the model
##
stepdrop1 <- function(object,  criterion="aic", steps=1000, k=2, alpha=0.05,
                      stat="wald",
                      headlong=FALSE, random=TRUE, details=1, trace=0,...){

  edlist <- getedges(object, complement=FALSE)
  edlist <- lapply(edlist, list)
  
  itcount <- 1
  
  t0 <- proc.time()
  repeat{
    ev <- evalInECC(object, edlist=edlist, headlong=headlong, random=random,
                      alpha=alpha,
                      stat=stat,
                      criterion=criterion, k=k, print=details-1)
    if (is.null(ev$optEdge) || itcount>steps)
      break()

    sss <- sprintf("Dropped edge: %30s   statistic %10.5f",
                   paste(unlist(ev$optEdge), collapse=" ~ "), ev$optStat)
    cat(sss, "\n")
    object  <- ev$optModel     


    edlist  <- getedges(object, complement=FALSE)
    edlist  <- lapply(edlist, list)
    itcount <- itcount+1

  }
  if (details>=1)
    cat("Time (stepwise forward): ", (proc.time()-t0)[3], "\n")
  
  return(object)
}


### PORTO - working version
stepjoin1 <- function(object, scope, type='ecc', criterion='aic', steps=1000, k=2,
                      alpha=0.05, stat="wald",
                      details=1,trace=0)
{
  criterion <- match.arg(criterion, c("test","aic","bic"))
  type      <- match.arg(type, c("ecc","vcc"))
  stat      <- match.arg(stat, c("wald","dev"))

  if (criterion %in% c("aic","bic")){
    alpha <- 0
  }
  
  if (details>=1){   
    typeStr <- ifelse (type=="vcc", "vertex colour classes", "edge colour classes")
    cat("Stepwise joining of", typeStr, "\n")
    cat(.stepDetailsString(criterion, stat, k, alpha),"\n")
  }   

  if (missing(scope)){
    scopecc <- getcc(object,type)
  } else{
    scopecc <- .addccnames(formula2names(scope),type=type)
  }

  currcc <- getcc(object,type)
  
  x<-sapply(scopecc, function(e){
    i <- match.containsLL2(e, currcc);
    i
  })
  names(x)<-NULL
  
  untouch     <- currcc[-x]
  len.untouch <- length(untouch)
  changelist  <- NULL
  modelChange <- FALSE

  stepcount   <- 1
  repeat {
    ##cat("Iteration - scopecc (at entry)\n");  print(scopecc)

    cctab <- comparecc(object, cc1=scopecc, cc2=scopecc, type=type, stat=stat,
                       details=details-1)

    tab   <- cctab$tab

    if (is.null(tab))
      return(object)
    
    statValue <- .extractStat(tab, criterion, k, direction="drop")
    idx       <- which.max(statValue)
    optstat   <- statValue[idx]
    if (details>=3)
      print(tab)

    if (optstat > alpha){
      modelChange <- TRUE
      ccid <- as.character(tab[idx,c("cc1","cc2")])
      cc   <- cctab$cc1[ccid]
      
      if (details>=1){
        cat("Joining", paste(ccid, collapse=' '),":\n",
            paste(formula2string(names2formula(cc)), collapse=";\n "),"\n")
        xxx <- tab[idx,-c(1:2)]
        cat(.printTestResult(xxx),"\n")
      }



      ## Delete those cc's which are to be joined
      scopecc[names(cc)] <- NULL;      

      if (length(scopecc)){ ## There are still cc's left
        scopecc <- .addccnames(.joincc(list(unlist(cc,recursive=FALSE)), scopecc), type)
      } else {
        scopecc <- .addccnames(list(unlist(cc, recursive=FALSE)),type)
      }
      
      ##cat("Updated scopecc:\n"); print(scopecc)

      if (length(untouch)>0){
        newcc <- .addccnames(unionL2L2(untouch, scopecc),type)
      } else {
        newcc <- scopecc
      }
      
      ##cat("CC for updated model:\n"); print(newcc)

      #KS <- fitInfo(object,"K")
      KS <- object$fitInfo$K
      ##cat("Updating model \n")
      switch(type,
             "ecc"={
               object  <- update(object, ecc=newcc, trace=trace, fit=FALSE)
             },
             "vcc"={
               object  <- update(object, vcc=newcc, trace=trace, fit=FALSE)
             })

      KS     <- findKinModel(object, KS, type=object$type,regularize=TRUE)      
      object <- fit(object, Kstart=KS)
      #object <- update(object, Kstart=KS, fit=TRUE)
      
      stepcount  <- stepcount + 1
      changelist <- c(changelist, cc)

      if (stepcount>steps)
        break()
    } else {
      if (details>=3){
        print(tab)
      }
      break()
    }
    
  }
  object$change <- changelist
  if (details>=1)
    cat("\n")
  if (modelChange)
    return(object)
}


      ##print(cc); print(cctab$cc1[names(cc)]); ##print(cctab$cc2[idx]); 
      
