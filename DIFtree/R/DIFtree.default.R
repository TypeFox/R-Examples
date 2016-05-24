#' @export

DIFtree.default <-
function(Y,
                            X,
                            model=c("Rasch","Logistic"),
                            type=c("udif","dif","nudif"),
                            alpha=0.05,
                            nperm=1000,
                            trace=FALSE,
                            penalize=FALSE,
                            ...){
  
  # check input 
  if(missing(Y)){
    stop("argument 'Y' is missing, with no default")
  }
  if(missing(X)){
    stop("argument 'X' is missing, with no default")
  }
  
  model <- match.arg(model)
  type  <- match.arg(type)
  
  # varify input 
  if(!all(as.vector(t(Y)) %in% c(0,1))){
    stop("Y must be a binary 0/1 matrix")
  }
  if(!is.data.frame(X)){
    stop("X must be of class 'data.frame'")
  }
  if(nrow(X)!=nrow(Y)){
    stop("Dimensions of X and Y don't match")
  }
  if(any(sapply(1:ncol(X),function(j) class(X[,j]))=="character")){
    stop("variable of class 'character' is not useful")
  }
  if(any(sapply(1:ncol(X),function(j) class(X[,j]))=="logical")){
    stop("variable of class 'logical' is not useful")
  }
  if(any(sapply(1:ncol(X), function(j) {
    is.numeric(X[,j]) && (mean(X[,j])<10e-6 | var(X[,j])==1)}))){
    stop("Don't use scaled covariates")
  }
  
  # set parameters 
  npersons <- nrow(Y)           # number of observations
  nitems   <- ncol(Y)           # number of items 
  y        <- as.vector(t(Y))   # response vector
  sumscore <- rowSums(Y)        # sum scores
  
  # modify DM_kov
  DM_kov <- prepareX(X)
  nvar   <- ncol(DM_kov)        # number of variables 
  
  # compute ordered values 
  ordered_values <- lapply(1:nvar, function(j){
    if(!all((DM_kov[,j] - round(DM_kov[,j])) == 0)){
      quantile(DM_kov[,j],seq(0.05,1,by=0.05))
    } else{
      unique(sort(DM_kov[,j]))
    }
  })
  names(ordered_values) <- names(DM_kov)
  
  n_levels <- sapply(1:nvar, function(j) length(ordered_values[[j]]))
  n_s      <- n_levels-1
  
  # call functions
  if(model=="Rasch"){
    output <- tree_Rasch(y,DM_kov,npersons,nitems,nvar,ordered_values,n_levels,n_s,alpha,nperm,trace,penalize)
    coefficients <- list("thetas"=output$thetas,
                         "betas_dif"=output$betas_dif,
                         "betas_nodif"=output$betas_nodif)
  } 
  if(model=="Logistic"){
    output <- tree_Logistic(Y,DM_kov,npersons,nitems,nvar,sumscore,ordered_values,n_levels,n_s,type,alpha,nperm,trace)
    if(type=="udif"){
      coefficients <- list("gammas_dif"=output$gammas_dif,
                           "gammas_nodif"=output$gammas_nodif,
                           "betas"=output$betas)
    } else{
      coefficients <- list("gammas_dif"=output$gammas_dif,
                           "gammas_nodif"=output$gammas_nodif,
                           "alphas_dif"=output$alphas_dif,
                           "alphas_nodif"=output$alphas_nodif)
    }
  }
  
  to_return <- list("splits"=output$splits,
                    "coefficients"=coefficients,
                    "pvalues"=output$pvalues,
                    "devs"=output$devs,
                    "crits"=output$crits,
                    "Y"=Y,
                    "X"=DM_kov,
                    "persons"=npersons,
                    "items"=nitems,
                    "call"=match.call())
                    
  class(to_return) <- "DIFtree"
  return(to_return)
    
}
