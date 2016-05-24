# =========== Parameter estimators ===========

#'@title Parameter estimation with separate regression analyses
#'  
#'@description Estimates the model parameters in \code{inner}, \code{reflective}, and
#'\code{formative} separately.
#'
#'@details \code{params.separate} estimates the statistical model described by \code{model}
#'
#'@template modelSpecification
#'  
#'@details
#'Model estimation proceeds as follows. The weights \code{W} and the
#'data covariance matrix \code{S} are used to calculate the composite covariance matrix \code{C}
#'and the indicator-composite covariance matrix \code{IC}. These are matrices are used to
#'separately estimate each of teh three model matrices \code{inner}, \code{reflective}, and
#'\code{formative}. This approach of estimating the parameter matrices separately is the
#'standard way of estimation in the PLS literature.
#'  
#'The default estimation approach is to estimate all parameters with a series of OLS 
#'regressions using \code{\link{estimator.ols}}.
#'  
#'@inheritParams matrixpls-common
#'  
#'@param parametersInner A function used to estimate the \code{inner} model matrix. The default is
#'  \code{\link{estimator.ols}}
#'  
#'@param parametersReflective A function used to estimate the \code{reflective} model matrix. The
#'  default is \code{\link{estimator.ols}}
#'  
#'@param parametersFormative A function used to estimate the \code{formative} model matrix. The
#'  default is \code{\link{estimator.ols}}
#'  
#'@param disattenuate \code{TRUE} or \code{FALSE} (default) indicating whether \code{C} should be
#'  disattenuated before applying \code{parametersInner}.
#'  
#'@param reliabilities A function that provides the reliability estimates for disattenuation. The
#'  default is \code{\link{reliability.weightLoadingProduct}}
#'  
#'@param ... All other arguments are passed through to \code{parametersInner},
#'\code{parametersReflective}, and\code{parametersFormative}
#'  
#'@return A named vector of parameter estimates.
#'  
#'@templateVar attributes params.separate,C IC inner reflective formative Q,c
#'@template attributes
#'  
#'@export

params.separate <- function(S, model, W, ...,
                              parametersInner = estimator.ols,
                              parametersReflective = estimator.ols,
                              parametersFormative = estimator.ols,
                              disattenuate = FALSE,
                              reliabilities = reliability.weightLoadingProduct){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  results <- c()
  
  # Calculate the composite covariance matrix
  C <- W %*% S %*% t(W)
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  reflectiveEstimates <- parametersReflective(S, nativeModel$reflective, W, ..., C = C, IC = t(IC))  
  
  if(disattenuate){
    
    Q <- reliabilities(S, reflectiveEstimates, W, ...)
    
    C <- C / sqrt(Q) %*% t(sqrt(Q))
    diag(C) <- 1
    
    # Fix the IC matrix. Start by replacing correlations with the corrected loadings
    tL <- t(reflectiveEstimates)
    IC[tL!=0] <- tL[tL!=0]
    
    # Disattenuate the remaining correlations
    IC[tL==0] <- (IC/sqrt(Q))[tL==0]
  }
  
  formativeEstimates <- parametersFormative(S, nativeModel$formative, W, ..., IC = IC)  
  innerEstimates <- parametersInner(S, nativeModel$inner, W, ..., C = C)  
  
  results <- c(estimatesMatrixToVector(innerEstimates, nativeModel$inner, "~"),
               estimatesMatrixToVector(reflectiveEstimates, nativeModel$reflective, "=~", reverse = TRUE),
               estimatesMatrixToVector(formativeEstimates, nativeModel$formative, "<~"))
  
  # Copy all non-standard attributes from the estimates objects
  
  for(object in list(innerEstimates, reflectiveEstimates, formativeEstimates)){
    for(a in setdiff(names(attributes(object)), c("dim", "dimnames", "class", "names"))){
      attr(results,a) <- attr(object,a)
      attr(W,a) <- NULL
    }
  }
  
  # Store these in the result object
  attr(results,"C") <- C
  attr(results,"IC") <- IC
  attr(results,"inner") <- innerEstimates
  attr(results,"reflective") <- reflectiveEstimates
  attr(results,"formative") <- formativeEstimates
  
  if(disattenuate){
    attr(results,"Q") <- Q
  }
  
  return(results)
}
