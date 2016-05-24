#' Compute the BIC of a given structure
#' @export
#' @param X the dataset
#' @param Z binary adjacency matrix of the structure (size p)
#' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
#' @param Bic_old BIC (vector) associated to Zold
#' @param Zold another structure with some common parts with Z (allows to compute only the differences, to be faster)
#' @param methode parameter for OLS (matrix inversion) methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
#' @param star boolean defining wether classical BIC or BIC* (over-penalized by a hierarchical uniform assumption to avoid over-learning)is computed
#' @return The vector of the BICs associated to each covariate (conditionnal distribution) according to the sub-regression structure.
#' 
#' @examples
#'    \dontrun{

#'    require(CorReg)
#' data=mixture_generator(n=15,p=5,valid=0)#dataset generation
#'Z=data$Z #binary adjacency matrix that describes correlations within the dataset
#'X=data$X_appr
#'Bic_null_vect=density_estimation(X=X)$BIC_vect
#'#Computes the BIC associated to each covariate (optional, BicZ can do it if not given as an input)
#'#computes the BIC associated to the structure
#'res=BicZ(X = X,Z = Z,Bic_null_vect=Bic_null_vect)
#'     }

BicZ<-function(X=X,Z=Z,Bic_null_vect=NULL,Bic_old=NULL,methode=1,Zold=NULL,star=FALSE){
   if(is.null(Bic_null_vect)){
      Bic_null_vect=density_estimation(X=X)$BIC_vect 
   }
   if(is.null(Zold)| is.null(Bic_old)){
      Zold=0*Z
      Bic_old=Bic_null_vect
   }
   res=.Call( "BicZ",as.matrix(X),Z,Bic_null_vect,Bic_old,methode,Zold, PACKAGE = "CorReg")
   if(star){
      res$BIC=sum(res$BIC)-ProbaZ(Z,star=TRUE)
   }
   return(res$BIC)
   
}