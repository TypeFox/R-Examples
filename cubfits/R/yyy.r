### These are functions to be exprted.

### Main
cubfits <- my.cubfits
cubpred <- my.cubpred
cubappr <- my.cubappr

### Utility
init.function <- my.init.function

estimatePhi <- function(fitlist, reu13.list, y.list, n.list,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  if(is.null(.cubfitsEnv$my.estimatePhiAll)){
    stop("init.function() is not called.")
  }
  .cubfitsEnv$my.estimatePhiAll(fitlist, reu13.list, y.list, n.list,
    E.Phi, lower.optim, upper.optim, lower.integrate, upper.integrate, control)
} # End of estimatePhi().

fitMultinom <- function(reu13.df, phi, y, n, phi.new = NULL, coefstart = NULL){
  if(is.null(.cubfitsEnv$my.fitMultinomAll)){
    stop("init.function() is not called.")
  }

  ret <- .cubfitsEnv$my.fitMultinomAll(reu13.df, phi, y, n, phi.new, coefstart)
  names(ret) <- names(reu13.df)

  tmp <- convert.b.to.bVec(ret)
  tmp <- convert.bVec.to.b(tmp, names(reu13.df), .cubfitsEnv$model)
  for(i.aa in 1:length(ret)){
    ret[[i.aa]]$coef.mat <- tmp[[i.aa]]$coef.mat
  }

  ret
} # End of fitMultinom().

