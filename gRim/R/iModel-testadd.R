#####################################################################
##
## Test addition of 'edge' to model 'object'
##
## If new model is decomposable and edge is in one clique only, then
## degrees of freedom are adjusted for sparsity
##
#####################################################################

testadd <- function(object, edge, k=2, details=1,...)
  UseMethod("testadd")

print.testadd <- function(x,  ...){

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

testadd.iModel <- function(object, edge, k=2, details=1, ...){

  edge <- rhsFormula2list(edge)[[1]]
  if (length(edge)!=2){
    stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))
  }

  model.type <- class(object)[1]

  if (is.null((amat <- list(...)$amat)))
    amat <- glist2adjMAT(object$glist)

  ## Is edge is in model? stop if not
  if (!subsetof(edge, colnames(amat)))
    stop(cat("variables:", edge, "not in model\n"))

  if (amat[edge[1],edge[2]]!=0)
    stop(cat("edge:", edge, "already in model\n"))

  ## Add edge to model
  ## FIXME: Fails if amat is sparse!
  amat[edge[1],edge[2]] <- amat[edge[2],edge[1]] <- 1L

  ## Is model graphical?
  cliq <- maxCliqueMAT(amat)$maxCliques
  isgraph <- length(cliq)==length(object$glist)

  ## Is model decomposable?
  isdecomp <- length(mcsMAT(amat))>0

  ## Is edge only in one clique in decomposable model?
  onlyinone <- FALSE
  if (isdecomp){
    idx   <- isin (cliq, edge, index=TRUE)
    onlyinone <- sum(idx) == 1
  }

  if (isdecomp && onlyinone && model.type %in% c("cModel","dModel")){
    ## If edge is in one clique only, do test in marginal table
    ##
    hostcq <- cliq[idx==1][[1]]
    set <- c(edge, setdiffPrim(hostcq, edge))

    switch(model.type,
           "cModel"={
             ans <- ciTest_mvn(list(cov=object$datainfo$S, n.obs=object$datainfo$n.obs),
                               set=set, ...)
           },
           "dModel"={
             ans <- ciTest_table(object$datainfo$data, set=set, ...)
           })
    extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
  } else {
    ## Make usual LR-test
    ##
    ob2   <- update(object, list(add.edge=edge))
    ans   <- .comparemodels(ob2,object)
    extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
  }
  extra2 <- list(aic=-(ans$statistic - k*ans$df), k=k)
  ret <- c(ans, extra, extra2)
  class(ret) <- "testadd"
  return(ret)
}
