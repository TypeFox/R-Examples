#' Derivatives of inverse of link function (internal use)
#' 
#' Computes derivatives of inverse of link functions (real estimates) with
#' respect to the beta parameters which are used for delta method computation
#' of the var-cov matrix of real parameters.
#' 
#' Note: that function was renamed to deriv_inverse.link to avoid S3 generic
#' class conflicts. The derivatives of the inverse of the link functions are
#' simple computations using the real values and the design matrix values.  The
#' body of the function is as follows:
#' 
#' \preformatted{switch(link, logit=x*real*(1-real), log=x*real,
#' loglog=-real*x*log(real), cloglog=-log(1-real)*x*(1-real), identity=x,
#' mlogit=x*real*(1-real)) }
#' 
#' @aliases deriv.inverse.link deriv_inverse.link
#' @param real Vector of values of real parameters
#' @param x Matrix of design values
#' @param link Type of link function (e.g., "logit")
#' @return Vector of derivative values computed at values of real parameters
#' @author Jeff Laake
#' @seealso \code{\link{inverse.link}}, \code{\link{compute.real}}
#' @keywords utility
deriv_inverse.link <-
function(real,x,link)
{
if(substr(link,1,6)=="mlogit" | substr(link,1,6)=="MLogit")link="MLogit"
real=as.vector(real)
switch(link,
logit=x*real*(1-real),
log=x*real,
loglog=-real*x*log(real),
cloglog=-log(1-real)*x*(1-real),
identity=x,
mlogit=x*real*(1-real),
sin=x*cos(asin(2*real-1))/2,
Sin=x*cos(asin(2*real-1))/2,
Logit=x*real*(1-real),
Log=x*real,
LogLog=-real*x*log(real),
CLogLog=-log(1-real)*x*(1-real),
Identity=x,
MLogit=x*real*(1-real)
)
}
