meta3 <- function(y, v, cluster, x, data, intercept.constraints=NULL, coef.constraints=NULL, 
                  RE2.constraints=NULL, RE2.lbound=1e-10,
                  RE3.constraints=NULL, RE3.lbound=1e-10,
                  intervals.type=c("z", "LB"), I2="I2q", R2=TRUE,
                  model.name="Meta analysis with ML",
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
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  ## Fixed factors with missing levels reported by David Stanley
  cluster <- eval(my.cluster, data, enclos = sys.frame(sys.parent()))
  if (!is.character(cluster))
    cluster <- as.character(cluster)
  ## check if there are missing data in cluster
  if (any(is.na(cluster)))
      stop("Missing values are not allowed in \"cluster\".\n")

  ## data in long format
  my.long <- data.frame(y, v, cluster)

  ## Add level-2 and level-3 predictors into the data set
  if (missing(x)) {
    no.x <- 0
    ## Indicator of either missing v or x
    miss.x <- is.na(v)
    } else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) no.x <- 1 else no.x <- ncol(x)   
    old.labels <- names(my.long)
    my.long <- data.frame(my.long, x)
    names(my.long) <- c(old.labels, paste("x", 1:no.x, sep=""))
    ## Indicator of either missing v or x
    miss.x <- apply(is.na(cbind(v,x)), 1, any)
  }
  
  ## Remove missing x. Missing y is fine as it will be handled by OpenMx.
  my.long <- my.long[!miss.x, ]
  ## Reorder data according to clusters
  my.long <- my.long[order(my.long$cluster), ]

  ## Convert long format to wide format as SEM uses wide format
  ## c() is required to convert matrix to vector when the data are balanced.
  ## c() is not required when the data are unbalanced.
  ## as.character() is required to null data with levels
  my.long$time <- c(unlist(sapply(split(my.long$y, as.character(my.long$cluster)), function(x) 1:length(x))))
  my.wide <- reshape(my.long, timevar="time", idvar=c("cluster"), direction="wide", sep="_")

  ## ## Replace "." with "_" since OpenMx does not allow "." as variable names
  ## names(my.wide) <- sub("\\.", "_", names(my.wide))

  ## maximum no. of data in level-2 unit
  k <- max(sapply( split(my.long$cluster, my.long$cluster), length))
  ## NA in v is due to NA in y in wide format
  ## Replace with 1e10 though it does not affect the analysis as NA in y
  temp <- my.wide[, paste("v", 1:k, sep="_")]
  temp[is.na(temp)] <- 1e10
  my.wide[, paste("v", 1:k, sep="_")] <- temp
  ## Missing indicator on y
  miss.y <- is.na(my.wide[, paste("y", 1:k, sep="_")])
  
  ## Prepare matrices
  if (is.null(intercept.constraints)) {
    Inter <- mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=0, labels="Intercept", name="Inter")
  } else {
    ## Convert intercept.constraints into a row matrix if it is not a matrix
    if (!is.matrix(intercept.constraints))
      intercept.constraints <- as.matrix(intercept.constraints)
    
    if (!all(dim(intercept.constraints)==c(1, 1)))
      stop("Dimensions of \"intercept.constraints\" are incorrect.\n")
    Inter <- as.mxMatrix(intercept.constraints, name="Inter")
  }
  ## Row vector of ones
  oneRow <- mxMatrix("Unit", nrow=1, ncol=k, name="oneRow")
  
  if ( no.x==0 ) {
    ## matrices with all zeroes; not actually used in the calculations
    Beta <- mxMatrix("Zero", nrow=1, ncol=1, name="Beta")
    mydata <- mxMatrix("Zero", nrow=1, ncol=k, name="mydata")
  } else {
    ## NA is not available for definition variable even y is NA.
    ## Replace NA with 0 in x. Since y is missing, 0 does not affect the results.
    for (i in 1:no.x) {
      temp <- my.wide[, paste(paste("x", i, sep=""), 1:k, sep="_")]
      temp[miss.y] <- 0
      my.wide[, paste(paste("x", i, sep=""), 1:k, sep="_")] <- temp
    }    
    mydata <- mxMatrix("Full", nrow=no.x, ncol=k, free=FALSE, name="mydata",
                      labels=c(outer(1:no.x, 1:k, function(x, y) paste("data.x", x,"_", y, sep = ""))))

    if (is.null(coef.constraints)) {
      Beta <- mxMatrix("Full", nrow=1, ncol=no.x, free=TRUE, values=0,
                       labels=paste("Slope_", 1:no.x, sep=""), name="Beta")
    } else {
    ## Convert coef.constraints into a row matrix if it is not a matrix
    if (!is.matrix(coef.constraints))
      coef.constraints <- t(as.matrix(coef.constraints))      
      Beta <- as.mxMatrix(coef.constraints, name="Beta")
    }    
  }
  
  if ( length(RE2.lbound) != 1 ) {
    warning("\"RE2.lbound\" should be a scalar.")
    RE2.lbound <- 1e-10
  }
  if ( length(RE3.lbound) != 1 ) {
    warning("\"RE3.lbound\" should be a scalar.")
    RE3.lbound <- 1e-10
  }

  if ( is.null(RE2.constraints) ) {
    ## Tau2 <- mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=0.01, labels="Tau2_2",                   
    ##                  lbound=RE2.lbound, name="Tau2")
    tau2 <- "0.01*Tau2_2"
  } else {
    ## Tau2 <- as.mxMatrix(RE2.constraints, name="Tau2", lbound=RE2.lbound)
    if (length(RE2.constraints)!=1) 
      stop("Length of \"RE2.constraints\" is not 1.\n")
    tau2 <- RE2.constraints
  }

  if ( is.null(RE3.constraints) ) {
    ## Tau3 <- mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=0.01, labels="Tau2_3",                   
    ##                  lbound=RE3.lbound, name="Tau3")
    tau3 <- "0.01*Tau2_3"
  } else {
    ## Tau3 <- as.mxMatrix(RE3.constraints, name="Tau3", lbound=RE3.lbound)
    if (length(RE3.constraints)!=1) 
      stop("Length of \"RE3.constraints\" is not 1.\n")
    tau3 <- RE3.constraints
  }

  Tau <- as.mxMatrix( Diag(c(tau2, tau3)), lbound=matrix(c(RE2.lbound,NA,NA,RE3.lbound), nrow=2, ncol=2), name="Tau")
  
  Id <- mxMatrix("Iden", nrow=k, ncol=k, name="Id")
  Ones <- mxMatrix("Unit", nrow=k, ncol=k, name="Ones") 
  ## conditional sampling variances
  V <- mxMatrix("Diag", nrow=k, ncol=k, free=FALSE, labels=paste("data.v", 1:k, sep="_"), name="V")
  ## expMean <- mxAlgebra( oneRow %x% inter + coeff2 %*% data2 + coeff3 %*% data3, name="expMean")
  expMean <- mxAlgebra( oneRow %x% Inter + Beta %*% mydata, name="expMean")
  Tau2 <- mxAlgebra(Tau[1,1], name="Tau2")
  Tau3 <- mxAlgebra(Tau[2,2], name="Tau3")
  expCov <- mxAlgebra( Ones %x% Tau3 + Id %x% Tau2 + V, name="expCov")

  ## Assume NA first
  mx0.fit <- NA
  if (no.x==0) {

    ## Calculate I2
    ## Based on Higgins and Thompson (2002), Eq. 9
    sum.w <- sum(1/my.long$v)
    sum.w2 <- sum(1/my.long$v^2)
    no.studies <- length(my.long$v)
    ## Typical V based on Q statistic
    qV <- matrix((no.studies-1)*sum.w/(sum.w^2-sum.w2), nrow=1, ncol=1)
    qV <- as.mxMatrix(qV)
    ## Typical V based on harmonic mean  
    hmV <- matrix(no.studies/sum.w, nrow=1, ncol=1)
    hmV <- as.mxMatrix(hmV)
    ## Typical V based on arithmatic mean
    amV <- matrix(mean(my.long$v))
    amV <- as.mxMatrix(amV)

    I2q_2 <- mxAlgebra( Tau2/(Tau2+Tau3+qV), name="I2q_2")
    I2q_3 <- mxAlgebra( Tau3/(Tau2+Tau3+qV), name="I2q_3")
    I2hm_2 <- mxAlgebra( Tau2/(Tau2+Tau3+hmV), name="I2hm_2")
    I2hm_3 <- mxAlgebra( Tau3/(Tau2+Tau3+hmV), name="I2hm_3")  
    I2am_2 <- mxAlgebra( Tau2/(Tau2+Tau3+amV), name="I2am_2")
    I2am_3 <- mxAlgebra( Tau3/(Tau2+Tau3+amV), name="I2am_3")
    
    ICC_2 <- mxAlgebra( Tau2/(Tau2+Tau3), name="ICC_2")
    ICC_3 <- mxAlgebra( Tau3/(Tau2+Tau3), name="ICC_3")
    I2 <- match.arg(I2, c("I2q", "I2hm", "I2am", "ICC"), several.ok=TRUE)

    ## I2 <- match.arg(I2, c("I2q", "I2hm", "I2am"), several.ok=TRUE)
    ci <- c(outer(I2, c("_2","_3"), paste, sep=""))

    ## Modified for OpenMx 2.0
    mx.model <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
                        Inter, Beta, mydata, Tau, Tau2, Tau3, V, expMean, expCov,
                        qV, hmV, amV, I2q_2, I2q_3, I2hm_2, I2hm_3, I2am_2, I2am_3, ICC_2, ICC_3,
                        mxExpectationNormal("expCov","expMean", dimnames=paste("y", 1:k, sep="_")),
                        mxFitFunctionML(),
                        mxCI(c("Inter","Beta","Tau", ci)))
  } else {
    ## no.x > 0
    ## Modified for OpenMx 2.0  
    mx.model <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
                       Inter, Beta, mydata, Tau, Tau2, Tau3, V, expMean, expCov,                        
                       mxExpectationNormal("expCov","expMean", dimnames=paste("y", 1:k, sep="_")),
                       mxFitFunctionML(),
                       mxCI(c("Inter","Beta","Tau")))

    ## Calculate R2
    if (R2) mx0.fit <- tryCatch( meta3(y=y, v=v, cluster=cluster, data=my.long, model.name="No predictor",
                                       suppressWarnings=TRUE, silent=TRUE), error = function(e) e )
   }

  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.model)
  
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = mx.fit <- tryCatch( mxRun(mx.model, intervals=FALSE, suppressWarnings=suppressWarnings,
                                  silent=silent, ...), error = function(e) e ),
    LB = mx.fit <- tryCatch( mxRun(mx.model, intervals=TRUE, suppressWarnings=suppressWarnings,
                                   silent=silent, ...), error = function(e) e ) )
 
  if (inherits(mx.fit, "error")) {
    cat("Error in running the mxModel:\n")
    stop(print(mx.fit))
  }

  ## my.long is complete data
  ## FIXME: remove miss.x (is it used by meta()?)
  out <- list(call=mf, I2=I2, R2=R2, data.wide=my.wide, data=my.long,
              no.y=1, no.x=no.x, miss.x=rep(FALSE, nrow(my.long)), mx.model=mx.model,
              mx.fit=mx.fit, mx0.fit=mx0.fit, intervals.type=intervals.type)
  class(out) <- c("meta", "meta3")
  return(out)
}

