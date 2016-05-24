#' Function to run a heuristic t test of the difference between two generalized correlations.
#' 
#' 
#' @param rxy {generalized correlation r*(x|y) where y is the kernel cause.}
#' @param ryx {generalized correlation r*(y|x) where x is the kernel cause.}
#' @param n {Sample size needed to determine the degrees of freedom for the t test.}
#' @importFrom psych paired.r
#' @return Prints the t statistics and p-values.  
#' @note This function requires Revele's R package called `psych' in memory. This test
#'   is known to be conservative (i.e., often fails to reject 
#'   the null hypothesis of zero difference between the two generalized 
#'   correlation coefficients.)
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @examples
#' 
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' g1=gmcxy_np(x,y)
#' n=length(x)
#' h1=heurist(g1$corxy,g1$coryx,n)
#' print(h1)
#' print(h1$t) #t statistic
#' print(h1$p) #p-value
#' @export

heurist=function(rxy,ryx,n){
tstat = paired.r(xy=rxy, xz=ryx, yz = min(rxy,ryx), n)
return(tstat)}
