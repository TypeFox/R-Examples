reml3 <- function(y, v, cluster, x, data, RE2.startvalue=0.1, RE2.lbound=1e-10,
                  RE3.startvalue=RE2.startvalue, RE3.lbound=RE2.lbound, RE.equal=FALSE,
                  intervals.type=c("z", "LB"), model.name="Variance component with REML",
                  suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  mf <- match.call()
  if (missing(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  my.y <- mf[[match("y", names(mf))]]
  my.v <- mf[[match("v", names(mf))]]
  my.cluster <- mf[[match("cluster", names(mf))]]  
  cluster <- as.character(eval(my.cluster, data, enclos = sys.frame(sys.parent())))
  ## check if there are missing data in cluster
  if (any(is.na(cluster)))
      stop("Missing values are not allowed in \"cluster\".\n")  
  cluster.order <- order(cluster)
  cluster <- cluster[cluster.order]
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  y <- y[cluster.order]
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  v <- v[cluster.order]
  no.studies <- length(y)

  ## Response matrix
  Y <- as.mxMatrix(matrix(y, ncol=1), name="Y")
  ## Intercept
  X <- matrix(1, ncol=1, nrow=no.studies)
  
  if (missing(x)) {
    no.x <- 0
  } else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) x <- x[cluster.order] else x <- x[cluster.order, ]
    X <- cbind(X, x)
    no.x <- ncol(X)-1
  }
  
  numStats <- no.studies-ncol(X)
  X <- as.mxMatrix(X)
  
  ## maximum no. of data in level-2 unit
  k <- lapply( split(cluster, cluster), length )

  ## Prepare V: fixed and known
  V <- as.mxMatrix( Diag(v), name="V" )

  ## Prepare Tau2 
  Tau_2 <- Diag(rep(paste(RE2.startvalue,"*Tau2_2", sep=""), no.studies))
  re2.lbound <- Diag(RE2.lbound, nrow=no.studies, ncol=no.studies)
  re2.lbound[re2.lbound==0] <- NA
  Tau_2 <- as.mxMatrix(Tau_2, name="Tau_2", lbound=re2.lbound)
  
  ## Prepare Tau3 
  Tau_3 <- lapply(k, function(x) {matrix(paste(RE3.startvalue,"*Tau2_3", sep=""), ncol=x, nrow=x)})
  Tau_3 <- bdiagMat(Tau_3)
  re3.lbound <- lapply(k, function(x) {matrix(RE3.lbound, ncol=x, nrow=x)})
  re3.lbound <- bdiagMat(re3.lbound)
  re3.lbound[re3.lbound==0] <- NA
  Tau_3 <- as.mxMatrix(Tau_3, name="Tau_3", lbound=re3.lbound)
 
  # Inverse of (V+Tau)
  W <- mxAlgebra(solve(V+Tau_2+Tau_3), name="W")
  alpha <- mxAlgebra( solve(t(X)%*%W%*%X) %*% t(X) %*% W %*% Y, name="alpha")
  # -2LL
  obj <- mxAlgebra( ( log(det(V+Tau_2+Tau_3)) + log(det(t(X)%*%W%*%X)) +
                       t(Y-X%*%alpha)%*%W%*%(Y-X%*%alpha) ), name="obj")

  # Creat model for REML
  ## reml.model <- mxModel(model=model.name, X, Y, V, Tau_2, Tau_3, W, alpha, obj,
  ##                       mxAlgebraObjective("obj"), mxCI(c("Tau2_2","Tau2_3")))

  # Constrain equal variances
  if (RE.equal) {
    reml.model <- mxModel(model=model.name, X, Y, V, Tau_2, Tau_3, W, alpha, obj,
                          mxFitFunctionAlgebra("obj"), mxCI(c("Tau2")))
    reml.model <- omxSetParameters(reml.model, labels=c("Tau2_2", "Tau2_3"), newlabels=c("Tau2", "Tau2"),
                                   values=c(RE2.startvalue, RE2.startvalue), lbound=c(RE2.lbound, RE2.lbound))
  } else {
    reml.model <- mxModel(model=model.name, X, Y, V, Tau_2, Tau_3, W, alpha, obj,
                          mxFitFunctionAlgebra(algebra="obj", numObs=no.studies, numStats=numStats), mxCI(c("Tau2_2","Tau2_3")))
  }

  ## ## Assuming NA first
  ## reml0.fit <- NA  
  ## if (no.x==0) {

  ##   ## Calculate I2
  ##   ## Based on Higgins and Thompson (2002), Eq. 9
  ##   sum.w <- sum(1/v)
  ##   sum.w2 <- sum(1/v^2)
  ##   ## Typical V based on Q statistic
  ##   qV <- matrix((no.studies-1)*sum.w/(sum.w^2-sum.w2), nrow=1, ncol=1)
  ##   qV <- as.mxMatrix(qV)
  ##   ## Typical V based on harmonic mean  
  ##   hmV <- matrix(no.studies/sum.w, nrow=1, ncol=1)
  ##   hmV <- as.mxMatrix(hmV)
  ##   ## Typical V based on arithmatic mean
  ##   amV <- matrix(mean(v))
  ##   amV <- as.mxMatrix(amV)

  ##   I2q_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]+qV), name="I2q_2")
  ##   I2q_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]+qV), name="I2q_3")
  ##   I2hm_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]+hmV), name="I2hm_2")
  ##   I2hm_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]+hmV), name="I2hm_3")  
  ##   I2am_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]+amV), name="I2am_2")
  ##   I2am_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]+amV), name="I2am_3") 
  ##   ICC_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]), name="ICC_2")
  ##   ICC_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]), name="ICC_3")

  ##   I2 <- match.arg(I2, c("I2q", "I2hm", "I2am", "ICC"), several.ok=TRUE)
  ##   ci <- c(outer(I2, c("_2","_3"), paste, sep=""))

  ##   reml.model <- mxModel(model=model.name, X, Y, V, Tau2, Tau3, W, alpha, obj,
  ##                         qV, hmV, amV, I2q_2, I2q_3, I2hm_2, I2hm_3, I2am_2, I2am_3, ICC_2, ICC_3,
  ##                         mxAlgebraObjective("obj"), mxCI(c("Tau_2","Tau3", ci)))
    
  ##   ## meta3 <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
  ##   ##                  inter, coeff, mydata, Tau2, Tau3, V, expMean, expCov,
  ##   ##                  qV, hmV, amV, I2q_2, I2q_3, I2hm_2, I2hm_3, I2am_2, I2am_3, ICC_2, ICC_3,
  ##   ##                  mxFIMLObjective("expCov","expMean", dimnames=paste("y_", 1:k, sep="")),
  ##   ##                  mxCI(c("inter","coeff","Tau2","Tau3", ci)))
  ## } else {
  ##   ## no.x > 0
  ##   reml.model <- mxModel(model=model.name, X, Y, V, Tau2, Tau3, W, alpha, obj,
  ##                         mxAlgebraObjective("obj"), mxCI(c("Tau_2","Tau3")))
  ##   ## meta3 <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
  ##   ##                  inter, coeff, mydata, Tau2, Tau3, V, expMean, expCov,
  ##   ##                  mxFIMLObjective("expCov","expMean", dimnames=paste("y_", 1:k, sep="")),
  ##   ##                  mxCI(c("inter","coeff","Tau2","Tau3")))

  ##   ## Calculate R2
  ##   if (R2) reml0.fit <- tryCatch( reml3(y=y, v=v, cluster=cluster, data=my.long, model.name="No predictor",
  ##                                  suppressWarnings=TRUE, silent=TRUE), error = function(e) e )    
  ## }

  ## Return mx model without running the analysis
  if (run==FALSE) return(reml.model)
  
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = mx.fit <- tryCatch( mxRun(reml.model, intervals=FALSE,
                                    suppressWarnings = suppressWarnings, silent=silent, ...), error = function(e) e ),
    LB = mx.fit <- tryCatch( mxRun(reml.model, intervals=TRUE,
                                     suppressWarnings = suppressWarnings, silent=silent, ...), error = function(e) e ) )
 
  if (inherits(mx.fit, "error")) {
    cat("Error in running the mxModel:\n")
    warning(print(mx.fit))
  }

  ## ## Ad-hoc: Add no. of studies and no. of observed statistics
  ## mx.fit@runstate$objectives[[1]]@numObs <- no.studies
  ## ## Ad-hoc: no. of observed statistics after removing the fixed-effects (p)
  ## mx.fit@runstate$objectives[[1]]@numStats <- numStats
  
  out <- list(call = mf, data=data, mx.model=reml.model, mx.fit=mx.fit, intervals.type=intervals.type, numObs=no.studies, numStats=numStats)
  class(out) <- c("reml", "reml3")
  return(out)
}
