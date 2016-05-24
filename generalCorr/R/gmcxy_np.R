#' Function to compute generalized correlation coefficients r*x|y and r*(y|x).
#'
#' This function uses the `np' package and assumes that there are no missing data.
#'
#'
#' @param x vector of x data
#' @param y vector of y data
#' @return
#' \item{corxy}{r*(x|y)  from regressing x on y, where y is the kernel cause.}
#' \item{coryx}{r*(y|x) from regressing y on x, where x is the cause.}
#' @note This is provided if the user want to avoid calling \code{kern}.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @references Vinod, H. D."Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics" in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#'
#' @references Vinod, H. D. "Matrix Algebra Topics in Statistics and Economics
#' Using R, Chapter 4 in Handbook of Statistics: Computational Statistics
#' with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @keywords kernel regression, asymmetric R*
#' @examples
#' \dontrun{
#' set.seed(34);x=sample(1:10);y=sample(2:11)
#' gmcxy_np(x,y)}
#'
#' ## The function is currently defined as
#' function (x, y)
#' {
#'     bw = npregbw(formula = x ~ y, tol = 0.1, ftol = 0.1)
#'     model = npreg(bws = bw, gradients = FALSE, residuals = TRUE)
#'     corxy = model$R2
#'     bw2 = npregbw(formula = y ~ x, tol = 0.1, ftol = 0.1)
#'     model2 = npreg(bws = bw2, gradients = FALSE, residuals = TRUE)
#'     coryx = model2$R2
#'     list(corxy = corxy, coryx = coryx)
#'   }
#' @export

gmcxy_np <-
function(x,y){
#np means we call the np library functions
bw=npregbw(formula=x~y,tol=0.1, ftol=0.1)
model=npreg(bws=bw, gradients=FALSE, residuals=TRUE)
corxy= model$R2
bw2=npregbw(formula=y~x,tol=0.1, ftol=0.1)
model2=npreg(bws=bw2, gradients=FALSE, residuals=TRUE)
coryx= model2$R2
list(corxy=corxy,coryx=coryx)}
