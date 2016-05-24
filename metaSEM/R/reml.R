reml <- function(y, v, x, data, RE.constraints=NULL, RE.startvalues=0.1, RE.lbound=1e-10,
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
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  
  if (is.vector(y)) {
    no.y <- 1
    no.studies <- length(y)
  } else {
    no.y <- ncol(y)
    no.studies <- nrow(y)
  }
  if (is.vector(v)) {
    no.v <- 1
  } else {
    no.v <- ncol(v)
  }
  if (missing(x)) {
    no.x <- 0
  } else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) {
      no.x <- 1
    } else {
      no.x <- ncol(x)
    }
  }
   
  # Check the dimensions of y and v
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
  v.labels <- vech(outer(1:no.y, 1:no.y, function(x, y) paste("v", x,"_", y, sep = "")))
  y.labels <- paste("y", 1:no.y, sep="")
  x.labels <- paste("x", 1:no.x, sep="")

  ## Create a design matrix: each row represents 1 study
  ## It is not used in the actural analysis.
  ## p: no. of predictors + intercept; remember to multiply by no.y
  if (missing(x)) {
    x.design <- matrix(1, ncol=1, nrow=no.studies)
    p <- no.y
    no.x <- 0
    input.df <- as.matrix(cbind(y,v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
  } else {
    x.design <- cbind(1, x)
    no.x <- ncol(x.design)
    p <- no.x*no.y
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, x.labels)) 
  }
  
  ## Create Y, a row vector of effect sizes
  ## nrow(Y)= no.y*no.studies
  Y <- c(t(y))
  ## A vector indicating missing values
  miss.vec <- is.na(Y)
  # miss.vec <- matrix( t(is.na(Y)), ncol=1 ) 
  # Remove missing values
  Y <- Y[!miss.vec]
  # Convert it into a column vector
  Y <- matrix(Y, ncol=1)
  # No. of total effect sizes after removing missing data
  no.es <- length(Y)
  # Ad-hoc: no. of observed statistics after removing the fixed-effects (p)
  numStats <- no.es-p
  Y <- as.mxMatrix(Y)
  
  # Function to create design matrix for a study, e.g., x_1=1,2,3, no.y=2
  # 1 0 2 0 3 0
  # 0 1 0 2 0 3
  fn1 <- function(x, no.y) {
    temp <- lapply(x, function(x, k){Diag(x=x, nrow=k, ncol=k)}, k=no.y)
    do.call(cbind, temp)
  }
  # temp: a list of design matrix per study
  temp <- lapply(split(x.design, 1:nrow(x.design)), fn1, no.y=no.y)
  # Convert the list into a design matrix
  # X: based on no. of effect sizes
  X <- do.call(rbind, temp)
  # A design matrix based on Y as a column vector
  X <- X[!miss.vec, , drop=FALSE ]
  X <- as.mxMatrix(X)
  
  # McCulloch (2003) Generalized linear mixed models. E-book p. 67
  # Searle, Casella, & McCulloch (1992). Variance Components. E-book p. 250
  # I - X (X'X)^-1 X'
  # crossprod(X)=X'X
  ## M <- diag(no.es) - X %&% solve( crossprod(X) )
  ## ##N <- diag(no.es) - X %*% solve( t(X)%*% X ) %*% t(X)
  ## M <- M[1:(no.es-p), ]  # remove redundant columns

  ## ## transformed effect size: A row vector as required
  ## y_star <- t( M %*% Y )  
  ## selVars <- paste("X", 1:(no.es-p), sep="")
  ## dimnames(y_star) <- list(NULL, selVars)
  ## M <- as.mxMatrix(M, name="M")
  
  ## matrix of lbound
  if (is.matrix(RE.lbound)) {
    
    if (!all(dim(RE.lbound)==c(no.y, no.y)))
      stop("Dimensions of \"RE.lbound\" are incorrect.")    
    } else {

      lbound <- matrix(NA, nrow=no.y, ncol=no.y)
      Diag(lbound) <- RE.lbound
    }  
  ## Convert it into a large matrix
  lbound <- bdiagRep(lbound, no.studies)
  lbound <- lbound[!miss.vec, !miss.vec]

  free <- bdiagRep( matrix(TRUE, nrow=no.y, ncol=no.y), no.studies )
  free <- free[!miss.vec, !miss.vec]
  ## free is a vector of logical values
  free <- as.logical(vech(free))
  
  ## Preparing the S matrix for covariance elements
  ## No predictor
  if (is.null(RE.constraints)) {
    
    if (is.matrix(RE.startvalues)) {
      # FIXME: test symmetry
      if (!all(dim(RE.startvalues)==c(no.y, no.y)))
        # stop instead of warning here
        stop("Dimensions of \"RE.startvalues\" are incorrect.")                    
      #values <- matrix(c(RE.startvalues), nrow=no.y, ncol=no.y)
    } else {      
      values <- Diag(x=RE.startvalues, nrow=no.y, ncol=no.y)
    }    
    # Large matrix
    values <- bdiagRep(values, no.studies)
    values <- values[!miss.vec, !miss.vec]
    
    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    # Large matrix
    Tau.labels <- bdiagRep(vec2symMat(Tau.labels), no.studies)
    Tau.labels <- Tau.labels[!miss.vec, !miss.vec]
    # Replace off diagonals with NA
    Tau.labels[Tau.labels=="0"] <- NA
    
    Tau <- mxMatrix("Symm", ncol=no.es, nrow=no.es, free=free, labels=vech(Tau.labels),
                    lbound=vech(lbound), values=vech(values), name="Tau")
  } else {
    ## Convert RE.constraints into a column matrix if it is not a matrix
    if (!is.matrix(RE.constraints))
      RE.constraints <- as.matrix(RE.constraints)
    
    if (!all(dim(RE.constraints)==c(no.y, no.y)))
      stop("Dimensions of \"RE.constraints\" are incorrect.")

    Tau <- bdiagRep(RE.constraints, no.studies)
    Tau <- Tau[!miss.vec, !miss.vec]
    Tau <- as.mxMatrix(Tau, lbound=lbound, name="Tau")
  }

  ## Known sampling variance matrix
  if (no.y==1) {
    V <- Diag(x=v, nrow=no.studies, ncol=no.studies)
  } else {
    V <- matrix2bdiag(v)
  }
  V <- V[!miss.vec, !miss.vec]
  V <- as.mxMatrix(V)

  # Inverse of (V+Tau)
  W <- mxAlgebra(solve(V+Tau), name="W")
  alpha <- mxAlgebra( solve(t(X)%*%W%*%X) %*% t(X) %*% W %*% Y, name="alpha")
  # -2LL
  obj <- mxAlgebra( ( log(det(V+Tau)) + log(det(t(X)%*%W%*%X)) +
                       t(Y-X%*%alpha)%*%W%*%(Y-X%*%alpha) ), name="obj")

  # Creat model for REML
  reml.model <- mxModel(model=model.name, X, Y, V, W, Tau, alpha, obj,
                   #mxData(observed=y_star, type="raw"),
                   mxFitFunctionAlgebra(algebra="obj", numObs=no.studies, numStats=numStats), mxCI("Tau"))

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

  ## Ad-hoc: Add no. of studies and no. of observed statistics
  ## mx.fit@runstate$objectives[[1]]@numObs <- no.studies
  ## mx.fit@runstate$objectives[[1]]@numStats <- numStats
  
  out <- list(call = mf, data=input.df, no.y=no.y, no.x=no.x, miss.vec=miss.vec, mx.model=reml.model,
              mx.fit=mx.fit, intervals.type=intervals.type, numObs=no.studies, numStats=numStats)
  class(out) <- "reml"
  return(out)
}
