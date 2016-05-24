#' @title edit PMIDs
#' @description add two sets of PMIDs together, or exclude one set PMIDs from another set of PMIDs.
#' @param x output of getPMIDs, or a set of PMIDs
#' @param y output of getPMIDs, or a set of PMIDs
#' @param method can be 'add' (default) or 'exclude'. see details.
#' @details when method is 'add', PMIDs in 'x' and 'y' will be combined. when method is 'exclude', PMIDs in 'y' will be excluded from 'x'. 
#' @seealso \code{\link{getPMIDs}}
#' @export
#' @examples
#' # pmid1=getPMIDs(author="Yan-Hui Fan",dFrom=2007,dTo=2013,n=10)
#' # rm1="22698742"
#' # pmids1=editPMIDs(x=pmid1,y=rm1,method="exclude")
#' 
#' # pmid2=getPMIDs(author="Yanhui Fan",dFrom=2007,dTo=2013,n=10)
#' # rm2="20576513"
#' # pmids2=editPMIDs(x=pmid2,y=rm2,method="exclude")
#' 
#' # pmids=editPMIDs(x=pmids1,y=pmids2,method="add")
editPMIDs <- function(x,y,method=c("add","exclude")){
  method=match.arg(method)
  if(any(sapply(x,str_length)!=8) || any(sapply(y,str_length)!=8)){
    stop("Each PMID should be 8 digits long")
  }
  pmids=c()
  if(method=="add"){
    pmids=c(x,y)
  }
  if(method=="exclude"){
    pmids=x[which(!x%in%y)]
  }
  pmids=unique(pmids)
  return(pmids)
}
