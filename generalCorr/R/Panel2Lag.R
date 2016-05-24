#' Function to compute a vector of 2 lagged values of a variable from panel data.
#'
#' The panel data have a set of time series for each entity (e.g. country)
#' arranged such that all time series data for one entitiy is together. The
#' data for the second entity should be below the entire data for first entity.
#' When a variable is lagged twice, special care is needed to insert NA's for
#' the first two time points (e.g. weeks) for each entity (country).
#'
#' @param ID {location of the column having time identities (e.g. the week number)}
#' @param xj {data on variable to be lagged linked to ID}
#' @return Vector containing  2 lagged values of xj.
#' @seealso A more general function \code{\link{PanelLag}} has examples.
#' @note This function is provided for convenient user modifications.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
## @examples
#'
#' @export

Panel2Lag <-
function(ID,xj){
nr=length(as.vector(ID))
outj=rep(NA,nr)
for (i in 1:nr){
if (ID[i]==1)  outj[i]=NA
if (ID[i]==2)  outj[i]=NA
if (ID[i]>2)  outj[i]=xj[i-2]}
return(outj)}
