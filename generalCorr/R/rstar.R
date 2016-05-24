#' Function to compute generalized correlation coefficients r*(x,y).
#'
#' Uses Vinod (2015) definition of generalized (asymmetric) correlation
#' coefficients.  It requires kernel regression of x on y obtained by using the `np' package.
#' It also reports usual Pearson correlation coefficient r and p-value for testing
#' the null hypothesis that (population r)=0.
#' @param x vector of data on the dependent variable
#' @param y data on the regressors which can be a matrix
#' @importFrom stats cor.test
#' @return Four items:
#' \item{corxy}{r*x|y or regressing x on y}
#' \item{coryx}{r*y|x or regressing y on x}
#' \item{pearson.r }{Pearson's product moment correlation coefficient}
#' \item{pv}{p-value for testing the Pearson r}
#' @note This function needs the kern function which in turn needs the np package.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also as \code{\link{gmcmtx0}}
#' @references Vinod, H. D."Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics" in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#'
#' @references Vinod, H. D. "Matrix Algebra Topics in Statistics and Economics
#' Using R", Chapter 4 in Handbook of Statistics: Computational Statistics
#' with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#'
#' @keywords asymmetric  p-values
#' @examples
#'
#' x=sample(1:30);y=sample(1:30); rstar(x,y)
#'
#' @export

rstar <-
function(x, y){
c1=cor.test(x,y)
sig=sign(c1$estimate)
pv=c1$p.value
pearson.r= c1$estimate
mod.1=kern(dep.y=x, reg.x=y)
mod.2=kern(dep.y=y, reg.x=x)
corxy= sqrt(mod.1$R2)*sig
coryx= sqrt(mod.2$R2)*sig
list(corxy=corxy, coryx=coryx,pearson.r=pearson.r,pv=pv)}
