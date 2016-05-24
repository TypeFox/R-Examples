#' @title Rasch Residual Factor Analysis
#' @export rfa
#' @description Calculation of the rasch residual factor analysis proposed by Wright (1996) and further discussed by Linacre (1998) to detect multidimensionality.
#'  
#' @details no details in the moment.
#' @param pers_obj an object of class \code{"pers"} as a result from function \code{\link{pers}}.
#' @param na_treat value to be assigned to residual cells which have missing data in the original response matrix. default is set to \code{na_treat=0} to set the residuals to 0, which implys that they are imputed as 'fitting data', i.e., zero residuals. This can attenuate contrasts (see. http://www.rasch.org/rmt/rmt142m.htm). An option is to set it to \code{na_treat=NA}.
#' @param tr a logical value indicating whether the data (the residual matrix) is transposed prior to  calculation. This would perform a person analysis rather than a item analysis. The default is set to item analysis.
#' @param use a character string as used in function \code{\link{cor}} or \code{\link{cov}}, giving a method for computing covariances or correlations in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". The default is set to \code{use="complete.obs"} which will exclude cases by listwise deletion to keep the correlation matrix positive definit.
#' @param res a character string defining which type of (rasch--) residual to analyze when computing covariances or correlations. This must be (exactly) one of the strings "sr" for score residuals , "stdr" for standardised residuals, "srsq" for score residuals squared, or "stdrsq" for standardised residuals squared. The default is set to \code{res="stdr"} refering to Linacre (1998).
#' @param method a character string as used in function \code{\link{cor}} or \code{\link{cov}}, indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman", can be abbreviated. The default is set to \code{method="pearson"}.
#' @param cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix.The default is set to \code{cor=TRUE} to use the correlation matrix.
#' @return An object of class \code{c("rfa","list")}.
#' @references Wright, B. D. (1996). Comparing Rasch measurement and factor analysis. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 3}(1), 3–24.
#' @references Linacre, J. M. (1998). Detecting multidimensionality: which residual data-type works best? \emph{Journal of outcome measurement, 2}, 266–283.
#' @exportClass rfa
#' @examples ######################
#' ########
#' data(bfiN) # loading reponse data
#' pers_obj <- pers(pair(bfiN))
#' result <- rfa(pers_obj)
#' summary(result)
#' plot(result)
#' #### 

rfa <- function(pers_obj, na_treat=0, tr=FALSE ,use="complete.obs", res="stdr", method="pearson", cor=TRUE){
  # returns a list with results from rasch residual factor (pc) analysis
  # func. by joerg-henrik heine jhheine(at)googlemail.com
  # needs func. expscore in i.expscore.R
  # see book by trevor bond & christine fox p. 253
  obj <- expscore(pers_obj, na_treat=na_treat) # calls internal function for residuals
  if (tr==TRUE) {obj <- lapply(obj,t)} # transpose all resid. matrices if tr == TRUE
  Eni <- obj$Eni # expected scores 
  Yni <- obj$Yni # "sr" - score residual  
  Zni <- obj$Zni # "stdr" - standardised residual 
  Y2ni <- obj$Y2ni # "srsq" - score Residual Squared
  Z2ni <- obj$Z2ni # "stdrsq" - standardised residual squared 
  # check of arguments 
  if( !(any(res==c("sr","stdr","srsq","stdrsq"))) ){stop("wrong type of residuals selected","\n", "check argument 'res'","\n")}
  if( class(cor)!="logical" ){stop("argument 'cor' must be of class \"logical\"","\n")}
  #assign the selected residual type
  if(res=="sr"){resi <- Yni}
  if(res=="stdr"){resi <- Zni}
  if(res=="srsq"){resi <- Y2ni}
  if(res=="stdrsq"){resi <- Z2ni}
  # compute either cor or cov matrix using other aruments 
  if(cor==TRUE){matr <- cor(x=resi, use=use, method=method)}
  if(cor==FALSE){matr <- cov(x=resi, use=use, method=method)}
  # make eigen value decomp.
  eig <- eigen(matr)
  # str(eig$values)
  rownames(eig$vectors) <- row.names(matr)
  ### component names
  Cname<-nchar(paste(dim(matr)[1]))
  colnames(eig$vectors) <-paste("Comp",formatC(1:dim(matr)[1], width = Cname, format = "d", flag = "0"),sep=".")
  
  eigen.values <- eig$values
  variance.total <- sum(eig$values) # variance to explain
  variance.proportion <- eig$values/sum(eig$values) # Proportion of Variance explained
  loadings <- eig$vectors # loadings
  
  erg <- list(pers_obj=pers_obj, pca=list(eigen.values=eigen.values, loadings=loadings, variance.proportion=variance.proportion, variance.total=variance.total), transposed=tr)
  class(erg) <- c("rfa","list")
  return(erg)
}
