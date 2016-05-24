#' Lotka's law estimation
#'
#' It estimates Lotka's law coefficients for scientific productivity (\cite{Lotka A.J., 1926})
#' @param results is an object of the class '\code{bibliometrix}' for which the analysis of the authors' dominance ranking is desired.
#' @return The function \code{lotka} returns a list of summary statistics of the Lotka's law estimation of an object of class \code{bibliometrix}.
#'
#' the list contains the following objects:
#' \tabular{lll}{
#' \code{Beta}  \tab   \tab Beta coefficient\cr
#' \code{C}   \tab   \tab Constant coefficient\cr
#' \code{R2} \tab   \tab Goodness of Fit\cr
#' \code{AuthorProd}    \tab   \tab Authors' Productivity frequency table}
#'
#' @examples
#' data(scientometrics)
#' results <- biblioAnalysis(scientometrics)
#' L=lotka(results)
#' L
#'
#' @seealso \code{\link{biblioAnalysis}} function for bibliometric analysis
#' @seealso \code{\link{summary}} method for class '\code{bibliometrix}'
#'


lotka<-function(results){

  # Author Productivity (LOTKA's LAW)
  Authors=results$Authors
  AUdf=data.frame(Authors)
  AuthorProd=aggregate(AUdf,by=list(AUdf$Freq),"length")
  AuthorProd[,2]=as.numeric(AuthorProd[,2])
  AuthorProd[,3]=AuthorProd[,2]/sum(AuthorProd[,2])
  names(AuthorProd)=c("N.Articles","N.Authors","Freq")
  LOTKA=(lm(log10(Freq)~log10(N.Articles),data=AuthorProd))

  L=list(Beta=abs(as.numeric(LOTKA$coeff[2])),C=as.numeric(LOTKA$coeff[1]), R2=summary(LOTKA)$r.squared,AuthorProd=AuthorProd)
  return(L)
}
