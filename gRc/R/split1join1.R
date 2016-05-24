
####
#### What do we do here???
####

split1 <- function(object, scope=NULL, type='ecc', details=1){
  logL0 <- logL(object)

  if (missing(scope)){
    ccl   <- getSlot(object, type)
    ccl   <- ccl[sapply(ccl, length)>1]
  } else {
    ccl   <- .addccnames(formula2names(scope), type=type)
  }
  
  if (length(ccl)==0)
    return(NULL)
  ans <- lapply(ccl, function(cc){
    switch(type, 
           "ecc"={
             mtmp <- update(object, splitecc=list(cc))
           },
           "vcc"={
             mtmp <- update(object, splitvcc=list(cc))
           }
         )
          
    list(cc2str(cc), 2*(logL(mtmp)-logL0), length(cc)-1)
  })
  
  ans <- .nestedList2df(ans)

  names(ans)    <- c("cc", "X2", "df")
  ccl           <- .addccnames(ccl,type)
  ans$cc        <- names(ccl)
  rownames(ans) <- 1:nrow(ans)
  ans           <- .addStat(ans,n=dataRep(object,"n"))

  attr(ans, "ccterms")<-ccl
  ans

  ans2 <- structure(list(tab=ans, cc=.addccnames(ccl,type), details=details),
                    class=c("statTable","data.frame"))
  ans2


}  


join1 <- function(object, scope=NULL, type='ecc',details=1, stat='wald'){

  if (missing(scope)){
    ccl   <- getSlot(object, type)
  } else {
    ccl   <- .addccnames(formula2names(scope), type=type)
  }

  ans <- comparecc(object, cc1=ccl, cc2=ccl, type=type, stat=stat, details=details)
  return(ans)
}

