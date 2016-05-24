print.mitml.summary <- function(x,...){
# print method for objects of class "summary.mitml"

  cl <- x$call
  vrs <- x$model
  itr <- x$iter
  ngr <- x$ngr
  mdr <- x$missing.rates
  conv <- x$conv

  # print general information
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")
  cat("\n")

  cat(formatC("Cluster variable:",width=-25), vrs$clus, sep=" ", collapse="\n")
  cat(formatC("Target variables:",width=-25), vrs$yvrs, collapse="\n")
  cat(formatC("Fixed effect predictors:",width=-25), vrs$pvrs, collapse="\n")
  cat(formatC("Random effect predictors:",width=-25), vrs$qvrs, collapse="\n")

  cat("\nPerformed", sprintf("%.0f",itr$burn), "burn-in iterations, and generated", sprintf("%.0f",itr$m),
      "imputed data sets,\neach", sprintf("%.0f",itr$iter), "iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", sprintf("%.0f",ngr), "groups.")},"\n")

  # print convergence diagnostics
  if(!is.null(conv)){

    for(cc in attr(conv,"stats")){

      cout <- matrix(c( sapply(conv, function(z) min(z[,cc])),
                sapply(conv, function(z) quantile(z[,cc],.25)),
                sapply(conv, function(z) mean(z[,cc])),
                sapply(conv, function(z) median(z[,cc])),
                sapply(conv, function(z) quantile(z[,cc],.75)),
                sapply(conv, function(z) max(z[,cc])) ), ncol=6 )
      rownames(cout) <- c("Beta:","Psi:","Sigma:")
      colnames(cout) <- c("Min","25%","Mean","Median","75%","Max")
      clab <- switch(cc, Rhat="\nPotential scale reduction (Rhat, imputation phase):\n",
                         SDprop="\nGoodness of approximation (imputation phase):\n")
      cat(clab,"\n")
      print.table(round(cout,3))

      clab <- switch(cc, Rhat="\nLargest potential scale reduction:\n",
                         SDprop="\nPoorest approximation:\n")
      cat(clab)
      maxb <- conv$beta[which.max(conv$beta[,cc]),]
      maxp <- conv$psi[which.max(conv$psi[,cc]),]
      maxs <- conv$sigma[which.max(conv$sigma[,cc]),]
      cat("Beta: [", paste(maxb[1:2],collapse=",") ,"], ",
          "Psi: [", paste(maxp[1:2],collapse=",") ,"], ",
          "Sigma: [", paste(maxs[1:2],collapse=",") ,"]\n", sep="")
    }
  }

  # missing data rates
  mdrout <- t(as.matrix(mdr))
  rownames(mdrout) <- "MD%"
  cat("\nMissing data per variable:\n")
  print.table(mdrout)

  cat("\n")

  invisible()
}

