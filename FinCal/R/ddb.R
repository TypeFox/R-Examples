#' Depreciation Expense Recognition -- double-declining balance (DDB), the most common declining balance method, which applies two times the straight-line rate to the declining balance.
#' 
#' @param cost cost of long-lived assets
#' @param rv   residual value of the long-lived assets at the end of its useful life. DDB does not explicitly use the asset's residual value in the calculations, but depreciation ends once the estimated residual value has been reached. If the asset is expected to have no residual value, the DB method will never fully depreciate it, so the DB method is typically changed to straight-line at some point in the asset's life.
#' @param t    length of the useful life
#' @seealso \code{\link{slde}}
#' @export
#' @examples
#' ddb(cost=1200,rv=200,t=5)
ddb <- function(cost,rv,t){
  if(t<2){
    stop("t should be larger than 1")
  }
  ddb=rep(0,t)
  ddb[1]=2*cost/t
  if(cost - ddb[1] <= rv){
    ddb[1]=cost-rv
  }else{
    cost=cost-ddb[1]
    for(i in 2:t){
      ddb[i]=2*cost/t
      if(cost - ddb[i] <= rv){
        ddb[i]=cost-rv
        break
      }else{
        cost = cost - ddb[i]
      }
    }
  }
  cbind(t=1:t,ddb=ddb)
}
