#' Inverse link functions (internal use)
#' 
#' Computes values of inverse of link functions for real estimates.
#' 
#' The inverse of the link function is the real parameter value. They are
#' simple functions of \code{X*Beta} where \code{X} is the design matrix values
#' and \code{Beta} is the vector of link function parameters. The body of the
#' function is as follows:
#' 
#' \preformatted{switch(link, logit=exp(x)/(1+exp(x)), log=exp(x),
#' loglog=exp(-exp(-x)), cloglog=1-exp(-exp(x)), identity=x,
#' mlogit=exp(x)/(1+sum(exp(x))) ) }
#' 
#' The \code{link="mlogit"} only works if the set of real parameters are
#' limited to those within the set of parameters with that specific link.  For
#' example, in POPAN, the \code{pent} parameters are of type "mlogit" so the
#' probabilities sum to 1.  However, if there are several groups then each
#' group will have a different set of \code{pent} parameters which are
#' identified by a different grouping of the "mlogit" parameters (i.e.,
#' "mlogit(1)" for group 1, "mlogit(2)" for group 2 etc).  Thus, in computing
#' real parameter values (see \code{\link{compute.real}}) which may have
#' varying links, those with "mlogit" are not used with this function using
#' \code{link="mlogit"}.  Instead, the link is temporarily altered to be of
#' type "log" (i.e., inverse=exp(x)) and then summed over sets with a common
#' value for "mlogit(j)" to construct the inverse for "mlogit" as
#' \code{exp(x)/(1+sum(exp(x))}.
#' 
#' @param x Matrix of design values multiplied by the vector of the beta
#' parameter values
#' @param link Type of link function (e.g., "logit")
#' @return Vector of real values computed from \code{x=X*Beta}
#' @author Jeff Laake
#' @seealso \code{\link{compute.real}},\code{\link{deriv_inverse.link}}
#' @keywords utility
inverse.link <-
function(x,link)
{
switch(link,
logit=1/(1+exp(-x)),
log=exp(x),
loglog=exp(-exp(-x)),
cloglog=1-exp(-exp(x)),
identity=x,
mlogit=1/(1+sum(exp(-x))),
probit=pnorm(x),
sin=(sin(x)+1)/2,
Logit=1/(1+exp(-x)),
Log=exp(x),
LogLog=exp(-exp(-x)),
CLogLog=1-exp(-exp(x)),
Identity=x,
MLogit=1/(1+sum(exp(-x))),
Sin=(sin(x)+1)/2
)
}
