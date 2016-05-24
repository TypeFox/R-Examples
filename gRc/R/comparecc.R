
## comparecc2 <- function(object, cc1=NULL, cc2=NULL, type='ecc', stat='wald',
##                        crit="aic",
##                        alpha=ifelse(crit=="aic", 0, 0.05),
##                        headlong = TRUE,
##                        details=1){


    
##   ## Relevant for Wald:
##   ##
##   if (stat=="wald"){
##     V <- vcov(object)
##     b <- coef(object)
##   }

##   ofs    <- ifelse(type=="ecc", length(object$vcc), 0)
##   cc      <- object[[type]]
##   #ccnames <- names(cc)
##   #lcc     <- length(cc)

  
##   if (is.null(cc1)){
##     cc1  <- cc;         ##cat("cc1 is NULL. Setting cc1 to all ccs\n")
##     cc1N <- object$intRep[[paste(type,"V",sep='')]]
##   }
##   if (is.null(cc2)){
##     cc2  <- cc;         ##cat("cc2 is NULL. Setting cc2 to all ccs\n")
##     cc2N <- object$intRep[[paste(type,"V",sep='')]]
##   }  
##   if (is.L(cc1)){
##     cc1 <- list(cc1);  ##cat("cc1 is a colour class. Changing to colour class list...\n")
##   }
##   if (is.L(cc2)){
##     cc2 <- list(cc2);  ##cat("cc2 is a colour class. Changing to colour class list...\n")
##   }

##   len1  <- length(cc1)
##   len2  <- length(cc2)
##   combs <- combn(len1, 2, simplify=FALSE) # FIXME: a hack...
  
##   optModel  <- object
##   optStat   <- -9999
##   optCC1    <- obtCC2 <- NULL
##   statlist  <- list()

##                                         #print(cc1N);  print(cc2N)
  
##   for (kk in 1:length(combs)){
##     ii <- combs[[kk]][1]
##     jj <- combs[[kk]][2]
##                                         #cat(sprintf("kk %d ii %d jj %d\n", kk, ii, jj))
##     if (!setequal(cc1N[[ii]], cc2N[[jj]])){
      
##       if (stat=="wald"){
##         cat("doing Wald")
##         cstat <- c(.findWald(ii+ofs, jj+ofs,  b,  V), 1)
##       } else{
##         cat("doing LR")
##         switch(type,
##                "ecc"={
##                  mtmp <- update(object, joinecc=list(cc1[[ii]], cc2[[jj]]))
##                },
##                "vcc"={
##                  mtmp <- update(object, joinvcc=list(cc1[[ii]], cc2[[jj]]))
##                })
##         cstat <- c(-2*(logL(mtmp)-logL(object)), 1) # Deviance for reduction
##       }
      
##       statval  <- compareModel(cstat, crit=crit)
##       cat("-------------\nstat:", statval,"\n")     
##       cat("cc1:", ii, paste(unlist(formula2string(names2formula(cc1[[ii]]))),collapse=" "),"\n")
##       cat("cc2:", jj, paste(unlist(formula2string(names2formula(cc2[[jj]]))),collapse=" "),"\n")

##       if (statval > max(optStat,alpha)){
##         optStat  <- stat
##         optii    <- ii
##         optjj    <- jj
##         optCC1   <- cc1[[ii]]
##         optCC2   <- cc1[[jj]]
##       }
##       statlist[[kk]] <- stat    
##       if (headlong && optStat>alpha)
##         break()
##     }
##   }
  
##   if (!is.null(optCC1)){
##     optModel <- update(object, joinecc=list(cc1[[ii]], cc2[[jj]]))
##   } else {
##     optModel <- object
##   }
  
##   ans <- list(optCC1=optCC1, optCC2=optCC2, optStat=optStat, optModel=optModel)

##   return(ans)
  
## }


comparecc <- function(object, cc1=NULL, cc2=NULL, type='ecc', stat='wald', details=1){
  cc      <- getSlot(object,type)
  ccnames <- names(cc)
  lcc     <- length(cc)

  if (details>=1)
    cat("Comparing colour classes of type:", type, "using statistic:", stat, "\n")
  
  if (lcc<=1)
    return(NULL)

  
  n   <- dataRep(object,"n")  
  type=match.arg(type, c("ecc", "vcc"))
  stat=match.arg(stat, c("wald", "dev"))
    
  ## Relevant for Wald:
  if (stat=="wald"){
    V <- vcov(object)
    b <- coef(object)
  }

  ofs   <- ifelse(type=="ecc", length(getSlot(object,"vcc")), 0)
  
  if (is.null(cc1)){
    cc1 <- cc;         ##cat("cc1 is NULL. Setting cc1 to all ccs\n")
  }
  if (is.null(cc2)){
    cc2 <- cc;         ##cat("cc2 is NULL. Setting cc2 to all ccs\n")
  }  
  if (is.L(cc1)){
    cc1 <- list(cc1);  ##cat("cc1 is a colour class. Changing to colour class list...\n")
  }
  if (is.L(cc2)){
    cc2 <- list(cc2);  ##cat("cc2 is a colour class. Changing to colour class list...\n")
  }

  cc1    <- .addccnames(.ccl2names(cc1, cc),type)
  cc2    <- .addccnames(.ccl2names(cc2, cc),type)
  lcc1   <- length(cc1)
  lcc2   <- length(cc2)
  cc1idx <- sapply(cc1, matchLL2, cc)
  cc2idx <- sapply(cc2, matchLL2, cc)
  lcc1   <- length(cc1)
  lcc2   <- length(cc2)
  ans    <- ccterms <- vector("list", lcc1*lcc2)
  ccnames1 <- names(.addccnames(cc1,type))
  ccnames2 <- names(.addccnames(cc2,type))

  ##<- cc1list <- cc2list
  rowidx  <- 0
  for (i in 1:lcc1){
    cc1visits  <- cc1idx[1:i]
    i1         <- cc1idx[i]
    i2         <- i1 + ofs

    #cci.name <- names(b)[i2]
    for (j in 1:lcc2){
      j1 <-cc2idx[j]
      if ((j1 %in% cc1visits) && (i1 %in% cc2idx)){
        ;    
      } else {
        j2 <- j1 + ofs
        #ccj.name  <- names(b)[j2]
        if (stat=="wald")
          w <- .findWald(i2, j2,  b,  V)
        else{
          switch(type,
                 "ecc"={
                   mtmp <- update(object, joinecc=cc[c(i1,j1)])
                 },
                 "vcc"={
                   mtmp <- update(object, joinvcc=cc[c(i1,j1)])
                 })
          w <- -2*(logL(mtmp)-logL(object)) # Deviance for reduction
        }
        rowidx <- rowidx+1

        ##print(list(i1,j1,cci.name, ccj.name, w))
        ##ans[[rowidx]] <- list(i1,j1,cci.name, ccj.name, w)
        #ans[[rowidx]] <- list(ccnames[i1],ccnames[j1], w)
        ans[[rowidx]] <- list(ccnames1[i],ccnames2[j], w)
        ccterms[[rowidx]] <- cc[c(i,j)]
      }
    }
  }

  ans        <- ans[!sapply(ans, is.null)]
  ccterms    <- ccterms[!sapply(ccterms, is.null)]
  ans        <- .nestedList2df(ans)  
  names(ans) <- c("cc1", "cc2", "X2")
  ans[,"cc1"] <- as.character(ans[,"cc1"])
  ans[,"cc2"] <- as.character(ans[,"cc2"])

  ans$df <- 1
  ans <- .addStat(ans, n=n)
  attr(ans, "ccterms") <- ccterms
  ans2 <- structure(list(tab=ans, cc1=.addccnames(cc1,type), cc2=.addccnames(cc2,type),
                         details=details),
                    class=c("statTable","data.frame"))

  ans2

}



  ##   ans$p  <- 1-pchisq(ans$X2, df=1)
  ##   ans$aic <- -ans$X2+2
  ##   ans$bic <- -ans$X2+log(n)*ans$df
