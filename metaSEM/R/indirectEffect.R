indirectEffect <- function(x, n, standardized=TRUE, direct.effect=TRUE, run=TRUE) {
  
  if (is.list(x)) {
    t(mapply(indirectEffect, x, n=n, MoreArgs=list(standardized=standardized, direct.effect=direct.effect)))
  } else {
    
  dimnames(x) <- list(c("y","m","x"),c("y","m","x"))

  if (standardized) {
    ## starting values for standard deviations
    stvalues <- sqrt(Diag(x))
    myMod <- mxModel("Standardized indirect effect", type="RAM", mxData(observed=x, type="cov", numObs=n),
                     manifestVars=c("y","m","x"),
                     ## Im: imaginary variables
                     latentVars=c("yLat","mLat","xLat","xIm1","mIm2","xIm2","varIm1","varIm2","covIm"),
                     mxPath(from="xLat", to="x", arrows=1, free=TRUE, values=stvalues[3], labels="xsd"),
                     mxPath(from="mLat", to="m", arrows=1, free=TRUE, values=stvalues[2], labels="msd"),                   
                     mxPath(from="yLat", to="y", arrows=1, free=TRUE, values=stvalues[1], labels="ysd"),
                     mxPath(from="xLat", arrows=2, free=FALSE, values=1),
                     mxPath(from="mLat", arrows=2, free=FALSE, values=0),
                     mxPath(from="yLat", arrows=2, free=FALSE, values=0),                   
                     mxPath(from="xIm1", arrows=2, free=FALSE, values=-1),
                     mxPath(from="xIm2", arrows=2, free=FALSE, values=-1),
                     mxPath(from="mIm2", arrows=2, free=FALSE, values=-1),
                     mxPath(from="varIm1", arrows=2, free=FALSE, values=1),
                     mxPath(from="varIm2", arrows=2, free=FALSE, values=1),
                     mxPath(from="covIm", arrows=2, free=FALSE, values=0),
                     mxPath(from="x", arrows=2, free=FALSE, values=0),
                     mxPath(from="m", arrows=2, free=FALSE, values=0),
                     mxPath(from="y", arrows=2, free=FALSE, values=0),
                     mxPath(from="xLat", to="mLat", arrows=1, free=TRUE, values=0.2, labels="a"),
                     mxPath(from="mLat", to="yLat", arrows=1, free=TRUE, values=0.2, labels="b"),
                     mxPath(from="xLat", to="yLat", arrows=1, free=TRUE, values=0.2, labels="c"),
                     mxPath(from=c("varIm1","xIm1"), to="mLat", free=c(FALSE, TRUE), values=c(1,0.2), labels=c(NA,"a")),
                     mxPath(from=c("varIm2","xIm2","mIm2"), to="yLat", free=c(FALSE, TRUE, TRUE),
                            values=c(1,0.2,0.2), labels=c(NA,"c","b")),
                     mxPath(from="xIm2", to="covIm", arrows=2, free=TRUE, values=0.2, labels="a"),
                     mxPath(from="covIm", to="mIm2", arrows=2, free=FALSE, values=-1))

  } else {
    myMod <- mxModel("Unstandardized indirect effect", type="RAM", mxData(observed=x, type="cov", numObs=n),
                     manifestVars=c("y","m","x"),
                     mxPath(from=c("x"), to=c("m"), arrows=1, free=TRUE, labels=c("a")),
                     mxPath(from=c("x","m"), to=c("y"), arrows=1, free=TRUE, labels=c("c","b")),
                     mxPath(from=c("x","m","y"), arrows=2, free=TRUE, values=c(1,0.2,0.2),
                            labels=c("xvar","mvar","yvar")))
  }
  ## No direct effect by setting c=0
  if (!direct.effect) myMod <- omxSetParameters(myMod, labels="c", free=FALSE, values=0) 

  ## Return mx model without running the analysis
  if (run==FALSE) return(myMod)
  
  my.fit <- mxRun(myMod, silent=TRUE, suppressWarnings=TRUE)
  ## my.parameters <- summary(my.fit)$parameters
  ## my.parameters$Estimate[my.parameters$name %in% "a"]
  a <- mxEval(a, my.fit) 
  b <- mxEval(b, my.fit) 
  indirect <- a*b

  # Fixed a bug that all elements have to be inverted before selecting some of them
  if (direct.effect) {
     acovS <- tryCatch( 2*solve(my.fit@output$calculatedHessian)[c("a","b","c"), c("a","b","c")], error = function(e) e)
     ## Differentiate of [a*b c]
     esDiff <- matrix(c(b, 0, a, 0, 0, 1), ncol=3)
  } else {
     acovS <- tryCatch( 2*solve(my.fit@output$calculatedHessian)[c("a","b"), c("a","b")], error = function(e) e)
     esDiff <- matrix(c(b, a), ncol=2)
  }
  
  if (inherits(acovS, "error")) {
    cat("Asymptotic covariance matrix of the estimates is not positive definite.\n")
    stop(print(acovS))
  }

  ## Acov of indirect (and direct effect)
  acovES <- esDiff %*% acovS %*% t(esDiff)

  if (direct.effect) {
     c(ind_eff=indirect, dir_eff=mxEval(c, my.fit), ind_var=acovES[1,1], ind_dir_cov=acovES[2,1],
       dir_var=acovES[2,2])
  } else {
     c(ind_eff=indirect, ind_var=acovES[1,1])
  }

  }  ## mapply
}

## es1: a*b;
## es2: c;
## esDiff: matrix([diff(es1, a), diff(es1, b), diff(es1, c)], 
##                [diff(es2, a), diff(es2, b), diff(es2, c)]);
## xVar: matrix ( [avar, abcov,accov], 
##                [abcov, bvar,bccov],
##                [accov, bccov, cvar]);
## esVar: esDiff . xVar . transpose(esDiff);
