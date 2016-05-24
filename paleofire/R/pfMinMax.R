#' MiniMax transformation of a charcoal serie
#' 
#' MiniMax transformation of a charcoal serie
#' 
#' 
#' @param serie Numeric, a vector of charcoal values.
#' @return \item{out}{A vector of minimax transformed values.}
#' @author O. Blarquez
#' @seealso \code{\link{pfTransform}}
#' @examples
#' 
#' ## Retrieve a site
#' ID=pfSiteSel(site_name=="Pas-de-Fond")
#' ## Or a group of sites (Western North America)
#' ID=pfSiteSel(id_region==c("WNA0"))
#' 
#' ## Extract data
#' A=pfExtract(ID)
#' 
#' ## Plot the first site raw charcoal data
#' par(mfrow=c(1,2))
#' plot(A[A[,1]==ID$id_site[1],3],A[A[,1]==ID$id_site[1],4],type="l",main=ID$site_name[1],
#'      xlab="Age",ylab="raw Char")
#' ## Minimax transformation
#' B=pfMinMax(A[A[,1]==ID$id_site[1],4])
#' ## Plot the first site Minimax transformed charcoal data
#' plot(A[A[,1]==ID$id_site[1],3],B,type="l",main=ID$site_name[1],
#'      xlab="Age",ylab="Minimax")
#' 
pfMinMax=function(serie){
  # Minimax transformation
serie <- (serie - min(serie,na.rm=TRUE))/(max(serie,na.rm=TRUE)-min(serie,na.rm=TRUE))
}
