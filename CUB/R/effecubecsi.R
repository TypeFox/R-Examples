# @title Auxiliary function for the log-likelihood estimation of CUBE models with covariates 
# only for the feeling component
# @description Compute the opposite of the scalar function that is maximized when running the 
# E-M algorithm for CUBE models with covariates only for the feeling component.
# @aliases effecubecsi
# @usage effecubecsi(param, ordinal, W, m)
# @param param Vector of length equal to NCOL(W) + 3 whose entries are the initial parameters estimates
# @param ordinal  Vector of ordinal responses
# @param W Matrix of the selected covariates for explaining the feeling component
# @param m Number of ordinal categories
#' @keywords internal 


effecubecsi <-
function(param,ordinal,W,m){
  q<-length(param)-3
  pai<-param[1]
  gama<-param[2:(q+2)]
  phi<-param[q+3]
  -loglikcubecsi(m,ordinal,W,pai,gama,phi)
}
