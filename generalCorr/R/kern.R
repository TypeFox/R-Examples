#' Function to run kernel regression with options for residuals and gradients.
#' 
#' 
#' @param dep.y {has the dependent (response) variable}
#' @param reg.x {has the regressor (stimulus) variable}
#' @param tol {=0.1 (default)}
#' @param ftol {=0.1 (default)}
#' @param gradients {set to TRUE if gradients computations are desired}
#' @param residuals {set to TRUE if residuals are desired}
#' @importFrom np npreg npregbw
#' @return Creates a model object `mod' 'containing kernel regression output.
#' Type names(mod) to find out the variety of outputs produced by `npreg' of `np' package.
#' @note This is a work horse for causal identification
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
## @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Vinod, H. D."Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics" in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048} 
#' @keywords amorphous partial derivative, apd
#' @keywords kernel regression residuals
#' @examples
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:50],ncol=2)
#' require(np)
#' k1=kern(x[,1],x[,2])
#' print(k1$R2) #prints the R square of the kernel regression
#' }
#' 
#' @export

kern <-
function(dep.y,reg.x,tol=0.1, ftol=0.1,
gradients=FALSE,residuals=FALSE){
gr=FALSE;resz=FALSE
if(gradients) gr=TRUE
if(residuals) resz=TRUE
#bandwidths for nonparametric regressions
bw=npregbw(ydat=as.vector(dep.y),xdat=reg.x,tol=tol, ftol=ftol)
mod=npreg(bws=bw, gradients=gr, residuals=resz)
return(mod)
}
