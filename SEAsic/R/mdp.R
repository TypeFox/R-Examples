#' Mean Difference for Pairs
#' @param x a column vector of scores on which the rsd is conditioned
#' @param g1 a column vector of equated scores based on a single subpopulation (aligned with elements in x)
#' @param g2 a column vector of equated scores based on a different single subpopulation (aligned with elements in x)
#' @param f a column vector of relative frequency associated with each raw score (can be based on either overall population or a subpopulation) (aligned with elements in x)
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized mean difference)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Kolen, M.J., & Brennan, R.L. (2004). Test equating, scaling, and linking: Methods and practices (2nd ed.). NY: Springer.}
#' }
#' @seealso \code{\link{dx}}
#' @examples
#' #Unstandardized MD for subpopulations 1 and 2 in the example data set, ex.data
#' mdp(x=ex.data[,1],g1=ex.data[,3],g2=ex.data[,4],f= ex.data[,8])
#' 
#' #Unstandardized MD for subpopulations 4 and 5 in the example data set, ex.data
#' mdp(x=ex.data[,1],g1=ex.data[,6],g2=ex.data[,7],f= ex.data[,8])
#' 
#' #Standardized MD for subpopulations 4 and 5 in the example data set, ex.data
#' mdp(x=ex.data[,1],g1=ex.data[,6],g2=ex.data[,7],f= ex.data[,8],s=4.2)
#' @note The equally weighted version of this index (Kolen & Brennan, 2004) can be obtained by inputting an f vector consisting of identical elements that sum to 1. For example, using f=c(rep(.047619,21)) with the example data set, ex.data. 
#' @return mean difference
#' @export

mdp<- function(x,g1,g2,f,s){
  if(missing(s))
    s <- 1
  
  md <- (sum((g1-g2)*f))/s  
  return(md)
}