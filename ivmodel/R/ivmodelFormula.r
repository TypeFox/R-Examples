ivmodelFormula <- function(formula,data,subset,na.action,
                            beta0=0,alpha=0.05,k=c(0,1), 
                            heteroSE = FALSE, clusterID = NULL, 
                            deltarange=NULL) {
  if(!inherits(formula,"formula")) {
  	stop("method is only for formula objects!")
  }
  # code gratefully stolen from ivreg() (package AER).
  if (missing(data)) 
    data <- environment(formula)
  mf = match.call()
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE      
  formula <- as.Formula(formula)

  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 
        1:2)
  has_dot <- function(formula) inherits(try(terms(formula),silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2)) {
      formula <- as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))
    }
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf, "numeric"); Y = matrix(as.numeric(Y),length(Y),1)
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf)
  
  mtZ <- delete.response(terms(formula, data = data, rhs = 2))
  Z <- model.matrix(mtZ, mf)
  
  if("(Intercept)" %in% colnames(X)) {
  	intercept=TRUE
  	X = X[,!(colnames(X) %in% "(Intercept)"),drop=FALSE]
  	Z = Z[,!(colnames(Z) %in% "(Intercept)"),drop=FALSE]
  	if(dim(Z)[2] < 1) stop("There aren't any instruments!")
  } else{
  	intercept=FALSE
  } 
  
  # Parse X and Z into D, X, and Z
  whichD = !(colnames(X) %in% colnames(Z))
  D = X[,whichD,drop=FALSE]
  if(dim(D)[2] != 1) {
  	stop("There must be exacty one endogenous variable! Please specify the formula correctly!")
  }
  unname(Z);unname(D); unname(Y)
  if(sum(!whichD) == 0) {
  	ivmodel(Y=Y,D=D,Z=Z,intercept=intercept,
  	                beta0=beta0,alpha=alpha,k=k,
  	                heteroSE=heteroSE,clusterID=clusterID,
  	                deltarange=deltarange)
  } else {
  	unname(X)
  	X = X[,!whichD,drop=FALSE]
  	whichZ = !(colnames(Z) %in% colnames(X))
    Z = Z[,whichZ,drop=FALSE]
  	ivmodel(Y=Y,D=D,Z=Z,intercept=intercept,
  	                beta0=beta0,alpha=alpha,k=k,
  	                heteroSE=heteroSE,clusterID=clusterID,
  	                deltarange=deltarange)
 }           	
}