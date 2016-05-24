#' Extract charcoal data for a list of sites
#' 
#' Extract charcoal data from an object returned by \code{\link{pfSiteSel}}
#' 
#' 
#' @param ID An object returned by \code{\link{pfSiteSel}}.
#' @return
#' 
#' \item{out}{A matrix of charcoal data with the following structure:
#' out[,1]=Site identifiers, out[,2]=Depths, out[,3]=Estimated ages,
#' out[,4]=Charcoal data.}
#' @author O. Blarquez
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
#' # Plot the first site raw charcoal data
#' plot(A[A[,1]==ID$id_site[1],3],A[A[,1]==ID$id_site[1],4],type="l",main=ID$site_name[1],
#'      xlab="Age",ylab="raw Char")
#' 
pfExtract=function(ID){
  ## Avoid no visible binding for global variable
  paleofiredata=NULL; rm(paleofiredata)
  
  # Extract data for sites
  data(paleofiredata, envir = environment())
  if(is.numeric(ID)){  Ext=paleofiredata[paleofiredata[,1] %in% ID,] } else Ext=paleofiredata[paleofiredata[,1] %in% ID$id_site,]
  return(Ext)
}
