################################################################
##
## Test deletion of 'edge' from model 'object'
##
## If model is decomposable and edge is in one clique only, then
## degrees of freedom are adjusted for sparsity
##
################################################################

testdelete <- function(object, edge, k=2, details=1,...)
  UseMethod("testdelete")

print.testdelete <- function(x,  ...){

  cat(sprintf("dev: %8.3f df:%3i p.value: %7.5f AIC(k=%3.1f): %6.1f edge: %s \n",
              x$statistic, x$df, x$p.value, x$k, x$aic, .toString(x$edge,':')))
  if (x$conmethod=="data.based"){
    if (x$details>0){
      cat("host: ", x$hostcq, "\n")
      cat("Notice: Test performed in saturated marginal model\n")
    }
  } else {
    if (x$details>0){
      cat("Notice: Test perfomed by comparing likelihood ratios\n")
    }
  }
  return(invisible(x))
}

testdelete.iModel <- function(object, edge, k=2, details=1, ...){

  ## cat("testdelete.iModel\n")
  edge <- rhsFormula2list(edge)[[1]]
  if (length(edge)!=2){
    stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))
  }
  
  model.type <- class(object)[1]
  ##cat(sprintf("testdelete.iModel model.type=%s\n", model.type))
  if (is.null((amat <- list(...)$amat)))
    amat <- glist2adjMAT(object$glist)
  
  ## Is edge is in model? stop if not
  if (!subsetof(edge, colnames(amat)))
    stop(cat("variables:", edge, "not in model\n"))
  if (amat[edge[1],edge[2]]!=1)
    stop(cat("edge:", edge, "not in model\n"))
  
  ## Is model graphical?
  cliq    <- maxCliqueMAT(amat)$maxCliques
  isgraph <- length(cliq)==length(object$glist)
  ## Is model decomposable?
  isdecomp <- length(mcsMAT(amat))>0
  
  ## Is edge only in one clique in decomposable model?
  onlyinone <- FALSE
  if (isdecomp){
    idx   <- isin (cliq, edge, index=TRUE)
    onlyinone <- sum(idx)==1
  }
  
  if (isdecomp && onlyinone && model.type %in% c("cModel","dModel")){
    ## If edge is in one clique only, do test in marginal table
    ##
    hostcq <- cliq[idx==1][[1]]
    set    <- c(edge, setdiffPrim(hostcq, edge))
    ##     cat(sprintf("CHK: edge: %s hostcq: %s\n", toString(edge), toString(hostcq)))

    switch(model.type,
           "cModel"={
             ans <- ciTest_mvn(list(cov=object$datainfo$S,n.obs=object$datainfo$n),
                               set=set, ...)
           },
           "dModel"={
             ans <- ciTest_table(object$datainfo$data, set=set, ...)
           })
    extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
  } else {
    ##cat(sprintf("CHK: Make usual LR-test - edge = {%s}\n", toString(edge)))
    ob2   <- update(object, list(drop.edge=edge))
    ans   <- .comparemodels(object,ob2)
    extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
  }
  extra2 <- list(aic=ans$statistic - k*ans$df, k=k)
  ret    <- c(ans, extra, extra2)
  class(ret) <- "testdelete"
  return(ret)
}

.comparemodels <- function(m1,m2) {
  devdiff <- 2*(m2$fitinfo$dev - m1$fitinfo$dev)
  dfdiff  <- m2$fitinfo$dimension['df'] - m1$fitinfo$dimension['df']
  #cat(sprintf(".comparemodels: lrtdiff=%f, dfdiff=%f\n", lrtdiff, dfdiff))
  ## NB: dev = K - 2logL
  list('statistic'=devdiff, 'df'=dfdiff, 'p.value'=1-pchisq(devdiff, dfdiff))
}




