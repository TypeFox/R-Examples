
## testIn:  Function which tests whether each edge in "edgeList" can be delete from model "object"
## testOut: Is similar but in the other direction.

## Author: Søren Højsgaard

## Input
## object   : imod-object
## edgeList : A list of edges; each edge is a vector

## Output
## A dataframe with test statistics (p-value or change in AIC), edges and logical
## telling if the edge can be deleted. 

## Known issues
## It is not tested whether edges in edgeList are in the model.
## Should check if edgeList is NULL


testEdges <- function(object, edgeMAT=NULL, ingraph=TRUE, criterion="aic", k=2, alpha=NULL,
                      headlong=FALSE, details=1,...){
  UseMethod("testEdges")
}

testEdges.iModel <- function(object, edgeMAT=NULL, ingraph=TRUE, criterion="aic", k=2, alpha=NULL,
                      headlong=FALSE, details=1,...){
  
  cl <- match.call()
  if (ingraph){
    cl[[1]] <- as.name("testInEdges")
  } else {
    cl[[1]] <- as.name("testOutEdges")
  }
  eval(cl)
}



testInEdges <- function(object, edgeMAT=NULL, criterion="aic", k=2, alpha=NULL, headlong=FALSE, details=1,...){

  criterion <- match.arg(criterion, c("aic","test"))
  
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

  if (is.null(alpha)){
    alpha <- if (criterion=="aic") 0 else 0.05
  }

  testFun <- if (headlong)
    .testInEdges_headlong
  else
    .testInEdges_all

  vn   <- object$varNames
  amat <- glist2adjMAT(object$glist, vn=vn)

  if (is.null(edgeMAT)){
    edgeMAT <- getInEdgesMAT(amat)
  }
  if (nrow(edgeMAT)==0)
    stop("There are no edges to test...\n")
  
  testFun(object, edgeMAT, comp.op=comp.op, crit.str=crit.str, alpha=alpha, k=k, amat=amat, vn=vn, ...)
}


.testInEdges_all <- function(object, edgeMAT, comp.op=`<`, crit.str="aic", alpha=0,k=2, amat, vn, ...)
{

  ##cat(".testInEdges_all\n")
  if (nrow(edgeMAT)==0)
    return(NULL)

  testMAT <- matrix(0,nrow=nrow(edgeMAT), ncol=4)
  colnames(testMAT) <- c("statistic","df","p.value","aic")
  indic <- rep.int(0,nrow(edgeMAT))

  for (ii in seq_len(nrow(edgeMAT))){
    ##print(edgeMAT[ii,])
    edgeTest <- testdelete(object, edgeMAT[ii,],k=k, amat=amat, ...)
    
    testMAT[ii,] <- as.numeric(edgeTest[c("statistic","df","p.value","aic")])
    curr.stat <- edgeTest[[crit.str]]

    if (comp.op(curr.stat, alpha)) {
      indic[ii] <- 1      
    }
    ##print("BBBBBBBBBBBBBBB")
  }

  ans <- cbind(
    as.data.frame(testMAT),
    as.data.frame(edgeMAT,stringsAsFactors=FALSE),
    action=c("-","+")[indic+1])
  return(ans)
}  

.testInEdges_headlong <- function(object, edgeMAT, comp.op=`<`, crit.str="aic", alpha=0,k=2, amat, vn, ...)
{
  if (nrow(edgeMAT)==0)
    return(NULL)

  testMAT <- matrix(0,nrow=nrow(edgeMAT), ncol=4)
  colnames(testMAT) <- c("statistic","df","p.value","aic")
  perm <- sample(nrow(edgeMAT))
  
  for (ii in seq_len(nrow(edgeMAT))){  
    edgeTest     <- testdelete(object, edgeMAT[perm[ii],],k=k, amat=amat, ...)
    testMAT[ii,] <- as.numeric(edgeTest[c("statistic","df","p.value","aic")])
    curr.stat    <- edgeTest[[crit.str]]
    if (comp.op( curr.stat, alpha)) {
      break      
    }
  }
  
  ans <- cbind(
               as.data.frame(testMAT[1:ii,,drop=FALSE]),
               as.data.frame(edgeMAT[perm[1:ii],,drop=FALSE],stringsAsFactors=FALSE),
               action=c(rep("-",ii-1),"+"))
  return(ans)
}  


testOutEdges <- function(object, edgeMAT=NULL, criterion="aic", k=2, alpha=NULL, headlong=FALSE, details=1,...){

  criterion <- match.arg(criterion, c("aic","test"))
  
  switch(criterion,
         "aic" ={opt.op    <- which.min
                 comp.op   <- `<`
                 outstring <- "change.AIC"
                 crit.str  <- "aic"},
         "test"={opt.op    <- which.max
                 comp.op   <- `<`
                 outstring <- "p.value"
                 crit.str  <- "p.value"
               })

  if (is.null(alpha)){
    alpha <- if (criterion=="aic") 0 else 0.05
  }

  testFun <- if (headlong)
    .testOutEdges_headlong
  else
    .testOutEdges_all

  vn   <- object$varNames
  amat <- glist2adjMAT(object$glist, vn=vn)

  if (is.null(edgeMAT)){
    edgeMAT <- getOutEdgesMAT(amat)
  }
  if (nrow(edgeMAT)==0)
    stop("There are no missing edges to test...\n")
  

  testFun(object, edgeMAT, comp.op=comp.op, crit.str=crit.str, alpha=alpha, k=k, amat=amat, vn=vn, ...)
}


.testOutEdges_all <- function(object, edgeMAT, comp.op=`<`, crit.str="aic", alpha=0,k=2, amat, ...)
{
  if (nrow(edgeMAT)==0)
    return(NULL)

  #cat(sprintf(" .testOutEdges_all k=%f\n", k))
  
  testMAT <- matrix(0,nrow=nrow(edgeMAT), ncol=4)
  colnames(testMAT) <- c("statistic","df","p.value","aic")
  indic <- rep.int(0,nrow(edgeMAT))
  
  for (ii in seq_len(nrow(edgeMAT))){
    edgeTest     <- testadd(object, edgeMAT[ii,],k=k, amat=amat,...)
    testMAT[ii,] <- as.numeric(edgeTest[c("statistic","df","p.value","aic")])
    curr.stat    <- edgeTest[[crit.str]]
    if (comp.op(curr.stat, alpha)) {
      indic[ii] <- 1      
    }
  }

  ans <- cbind(
    as.data.frame(testMAT),
    as.data.frame(edgeMAT,stringsAsFactors=FALSE),
    action=c("-","+")[indic+1])
  return(ans)
}  


.testOutEdges_headlong <- function(object, edgeMAT, comp.op=`<`, crit.str="aic", alpha=0,k=2, amat, vn, ...)
{
  if (nrow(edgeMAT)==0)
    return(NULL)

  testMAT <- matrix(0,nrow=nrow(edgeMAT), ncol=4)
  colnames(testMAT) <- c("statistic","df","p.value","aic")
  
  perm <- sample(nrow(edgeMAT))
  for (ii in seq_len(nrow(edgeMAT))){  
    edgeTest <- testadd(object, edgeMAT[perm[ii],],k=k, amat=amat, ...)
    testMAT[ii,] <- as.numeric(edgeTest[c("statistic","df","p.value","aic")])
    curr.stat <- edgeTest[[crit.str]]
    if (comp.op( curr.stat, alpha)) {
      break      
    }
  }
  
  ans <- cbind(
               as.data.frame(testMAT[1:ii,,drop=FALSE]),
               as.data.frame(edgeMAT[perm[1:ii],,drop=FALSE],stringsAsFactors=FALSE),
               action=c(rep("-",ii-1),"+"))
  return(ans)
}  

