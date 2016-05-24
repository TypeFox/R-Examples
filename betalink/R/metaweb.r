#' @title Returns a metaweb given a list of networks
#' @description
#' Given a list of networks, this function returns the metaweb
#'
#' @param n a \code{list} of graphs
#' @export
metaweb <- function(n){
   if(!is.list(n))
   {
      warning("n must be given as a list\nI have converted it for you\nCHECK THE RESULTS")
      n <- list(n)
   }
   n <- name_networks(n)
   if(length(n)==1)
   {
      warning("There is a single network in the data you passed")
      return(n)
   }
   M <- n[[1]]
   for(i in c(2:length(n))) M <- M + n[[i]]
   return(M)
}
