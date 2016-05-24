#' Calculates error, long term power and pass/fail criteria for CI obtained from any method
#' @param n - Number of trials
#' @param LL - Lower limit
#' @param UL - Upper limit
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of intervals obtained from any method using error due to
#' the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family General methods for error
#' @examples
#' LL=c(0,0.01,0.0734,0.18237,0.3344,0.5492)		#Lower and Upper Limits
#' UL=c(0.4507,0.6655,0.8176,0.9265,0.9899,1)
#' n= 5; alp=0.05;phi=0.05; f=-2
#' errGEN(n,LL,UL,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errGEN<-function(n,LL,UL,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(LL)) stop("'Lower limit' is missing")
  if (missing(UL)) stop("'Upper Limit' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  if ((class(LL) != "integer") & (class(LL) != "numeric") || any(LL < 0)) stop("'LL' has to be a set of positive numeric vectors")
  if ((class(UL) != "integer") & (class(UL) != "numeric") || any(UL < 0)) stop("'UL' has to be a set of positive numeric vectors")
  if (length(LL) <= n ) stop("Length of vector LL has to be greater than n")
  if (length(UL) <= n ) stop("Length of vector UL has to be greater than n")
  if (any(LL[0:n+1] > UL[0:n+1] )) stop("LL value have to be lower than the corrosponding UL value")

x=0:n
k=n+1
####INITIALIZATION
delalp=0
alpstar=0
thetactr=0

for(m in 1:k)
{
if(phi > UL[m] || phi<LL[m])
{
thetactr=thetactr+1
alpstar[m]=dbinom(x[m],n,phi)
} else alpstar[m] = 0
}
delalp=round((alp-sum(alpstar))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalp < f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalp,theta,Fail_Pass))
}
