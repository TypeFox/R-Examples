#' calculate weighted average shares -- weighted average number of common shares
#'
#' @param ns n x 1 vector vector of number of shares
#' @param nm n x 1 vector vector of number of months relate to ns
#' @seealso \code{\link{EPS}}
#' @seealso \code{\link{diluted.EPS}}
#' @export
#' @examples
#' s=c(10000,2000);m=c(12,6);was(ns=s,nm=m)
#' 
#' s=c(11000,4400,-3000);m=c(12,9,4);was(ns=s,nm=m)
was <- function(ns,nm){
  m=length(ns)
  n=length(nm)
  sum=0
  if(m==n){
    for(i in 1:m){
      sum=sum+ns[i]*nm[i]
    }
  }else{
    stop("length of ns and nm must be equal")
  }
  sum=sum/12
  return(sum)
}
