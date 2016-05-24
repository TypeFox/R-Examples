#' Function for computing a vector of one-lagged values of xj, a variable from panel data.
#'
#' Panel data have a set of time series for each entity (e.g. country)
#' arranged such that all time series data for one entitiy is together, and the
#' data for the second entity should be below the entire data for first entity
#' and so on for entities. In such a data setup,
#' When a variable is lagged once, special care is needed to insert an NA for
#' the first time point in the data (e.g. week) for each entity.
#'
#' @param ID {location of the column having time identities (e.g. week number).}
#' @param xj {data on variable to be lagged linked with the ID.}
#' @param lag {Number of lags desired (lag=1 is the default).}
#' @return Vector containing one-lagged values of variable xj.
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @examples
#'
#' \dontrun{
#' indiv=gl(6,12,labels=LETTERS[1:6])  
#' #creates A,A,A 12 times B B B also 12 times etc.
#' set.seed(99);cost=sample(30:90, 72, replace=TRUE)
#' revenu=sample(50:110, 72, replace=TRUE); month=rep(1:12,6)
#' df=data.frame(indiv,month,cost,revenu);head(df);tail(df)
#' L2cost=PanelLag(ID=month,xj=df[,"cost"], lag=2)
#' head(L2cost)
#' tail(L2cost)
#' 
#' gmcmtx0(cbind(revenu,cost,L2cost))
#' 
#' gmcxy_np(revenu,cost)
#' }
#'
#' @export

PanelLag <-
function(ID,xj,lag=1){
nr=length(as.vector(ID))
outj=rep(NA,nr)
ID=as.numeric(ID)
for (i in 1:nr){
if (ID[i]<=lag)  outj[i]=NA
if (ID[i]>lag)  outj[i]=xj[i-lag]}
return(outj)}
