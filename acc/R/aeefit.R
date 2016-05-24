#' @export
#' @importFrom utils head tail 
#' @importFrom plyr ddply
#' @importFrom nleqslv nleqslv
#' @importFrom stats model.matrix
#' @importFrom methods getClass



##############################################################################
# User's Main Function
##############################################################################
aeefit <- function(formula, data, weight=NULL, control=list()) {
  
  method <- "AEE"
  se <- "Sandwich"
  Call <- match.call()
  fm <- formula
  
  if (is.null(weight)) {
    weight <- rep(1,nrow(data[!duplicated(data$ID),]))
  } 
  
  # A PanelSurv object
  obj <- eval(formula[[2]], data)
  
  # Combine respones data frame and covariate data frame (remove intercept column)
  # Multiple rows per subject
  formula[[2]] <- NULL
  
  if (formula == ~ 1) {
    DF <- cbind(obj$psDF, zero=0)
  } else {
    DF <- cbind(obj$psDF, model.matrix(formula, data))[, -4]
  }
  
  DF <- DF[order(DF$ID, DF$time), ]
  
  # Design matrix, one row per subject
  X <- as.matrix(ddply(DF, "ID", head, n=1)[, -c(1:3)])
  
  # Create an Engine object
  engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
  engine <- do.call("new", c(list(Class=method), engine.control))
  
  if (length(engine@betaInit) == 1 & ncol(X) > 1)
    engine@betaInit <- rep(engine@betaInit, ncol(X))
  if (length(engine@betaInit) > 1 & length(engine@betaInit) != ncol(X))
    stop("Invalid length of initial beta values!")
  
  # Create a StdErr object
  #if(se == "NULL"){
  #  stdErr <- NULL}
  #if(se != "NULL"){
    stdErr.control <- control[names(control) %in% names(attr(getClass(se), "slots"))]
    stdErr <- do.call("new", c(list(Class=se), stdErr.control))
  #}
  
  
  fit <- doPanelFit.AEE.Sandwich(DF=DF, panelMatrix=obj$panelMatrix, timeGrid=obj$timeGrid,
                                 X=X, engine=engine,weight=weight)
  
  ret = list(formula=fm, beta=fit$beta, 
        baseline=fit$baseline,
        timeGrid=fit$timeGrid,
        lambda=fit$lambda,
        convergence=fit$convergence,
        iter=fit$iter,
        betaSE=fit$betaSE,
        betaVar=fit$betaVar,
        baselineSE=fit$baselineSE)
  
  class(ret) <- "aeefit"
  ret
  
}

