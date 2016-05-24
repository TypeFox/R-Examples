#' @title CovLP
#'
#' @export
#'
setClass("CovDepthWeighted", representation(depth = "character"), contains="CovRobust")

#'@title CovLp
#'
#'@description Weighted by \eqn{L^p} depth (outlyingness) multivariate location and scatter estimators.
#'
#'  @param x The data as a matrix or data frame. If it is a matrix or data frame, then each row is viewed as one multivariate observation.
#'  @param pdim The parameter of the weighted \eqn{L^pdim} depth
#'  @param la parameter of a simple weight function w=a*x+b
#'  @param lb parameter of a simple weight function w=a*x+b
#'
#'  
#'  @return loc: Robust Estimate of Location:
#'  @return cov:  Robust Estimate of Covariance:
#'  @return Returns depth weighted covariance matrix.
#'  
#'
#'  @details
#'  
#'  Using depth function one can define a depth-weighted location and scatter estimators. In case of location estimator we have                                                                                                                                                                                                                                                                                      \deqn{ L(F)={\int{{x}{{w}_{1}}(D({x},F))dF({x})}}/{{{w}_{1}}(D({x},F))dF({x})}, }                                                                                                                                                                                                                                                                                        Subsequently, a depth-weighted scatter estimator is defined as                                                                                                                                                                                                                                                                                      \deqn{ S(F)=\frac{\int{({x}-L(F)){{({x}-L(F))}^{T}}{{w}_{2}}(D({x},F))dF({x})}}{\int{{{w}_{2}}(D({x},F))dF({x})}}, }                                                                                                                                                                                                                                                                                        where  \eqn{ {{w}_{2}}(\cdot ) }  is a suitable weight function that can be different from  \eqn{ {{w}_{1}}(\cdot ) } .                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
#'  
#'  The \pkg{DepthProc} package offers these estimators for weighted  \eqn{ {L}^{p} } depth. Note that \eqn{ L(\cdot ) }  and  \eqn{ S(\cdot ) }  include multivariate versions of trimmed means and covariance matrices. Their sample counterparts take the form                                                                                                                                                                                                                                                                                      \deqn{ {{T}_{WD}}({{{X}}^{n}})={\sum\limits_{i=1}^{n}{{{d}_{i}}{{X}_{i}}}}/{\sum\limits_{i=1}^{n}{{{d}_{i}}}} , }                                                                                                                                                                                                                                                                                         \deqn{ DIS({{{X}}^{n}})=\frac{\sum\limits_{i=1}^{n}{{{d}_{i}}\left( {{{X}}_{i}}-{{T}_{WD}}({{{X}}^{n}}) \right){{\left( {{{X}}_{i}}-{{T}_{WD}}({{{X}}^{n}})                                                     \right)}^{T}}}}{\sum\limits_{i=1}^{n}{{{d}_{i}}}}, }    where  \eqn{ {{d}_{i}} }  are sample depth weights,   \eqn{ {{w}_{1}}(x)={{w}_{2}}(x)=x } .
#'  
#'  @author Daniel Kosiorowski and Zygmunt Zawadzki from Cracow University of Economics.
#'  
#'  @export
#'  @seealso \code{\link{depthContour}} and \code{\link{depthPersp}} for depth graphics.
#'  
#'  @examples
#'  x = mvrnorm(n = 100, mu = c(0,0), Sigma = 3*diag(2))
#'  cov_x = CovLP(x, 2, 1, 1)
#'  
#'  # EXAMPLE 2
#'  data(under5.mort,inf.mort,maesles.imm)
#'  data1990 = na.omit(cbind(under5.mort[,1],inf.mort[,1],maesles.imm[,1]))
#'  CovLP(data1990)
#'  
#'  
#'  @keywords
#'  multivariate
#'  nonparametric
#'  robust
#'  depth function
#'  
CovLP = function(x, pdim=2, la=1, lb=1)
{
  if(is.data.frame(x)) x = as.matrix(x)
  cov = CovLPCPP(x, pdim, la, lb)
  center = depthMedian(x, method = "LP", pdim = pdim, la = la, lb = lb)
  
  method = "Depth Weighted Estimator"
  new("CovDepthWeighted", cov = cov,
      center = center,
      det = det(cov),
      n.obs = nrow(x),
      X = x,
      method = method,
      call = match.call())
}


