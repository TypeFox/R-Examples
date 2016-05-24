

.depthLocal = function(u, X, beta, depth1, depth2, ...)
{
  ncol = ncol(X)
  nrow = nrow(X)
  
  ### Fix for dim 1
  tmp = t(apply(X, 1, function(k) 2*u-k))
  if(ncol(tmp) != ncol(X)) tmp = t(tmp)
  
  symDATA = rbind(X, tmp)
  depths = as.numeric(depth(X, symDATA, method=depth1, ...))
  quan = quantile(depths, probs=1-beta)
  Rset = as.matrix(X[ signif(depths,digits=6) >= signif(quan,digits=6), ])
  as.numeric(depth(u, Rset, method=depth2, ...))
}


#'@title Local depth
#'
#'@description Computes local version of depth according to proposals of Paindaveine and Van Bever - see referencess.
#'
#'@export
#'
#'  @param u Numerical vector or matrix whose depth is to be calculated. Dimension has to be the same as that of the observations.
#'  @param X The data as a matrix, data frame. If it is a matrix or data frame, then each row is viewed as one multivariate observation.
#'  @param beta cutoff value for neighbourhood
#'  @param depth1 depth method for symmetrised data
#'  @param depth2 depth method for calculation depth of given point
#'  @param name name for this data set - it will be used on plots.
#'  @param ... additional parameters passed to depth1 and depth2
#'  
#'  
#'  @details
#'  
#'  A successful concept of {local depth } was proposed by Paidaveine and Van Bever (2012) . For defining {a neighbourhood} of a point authors proposed using idea of {symmetrisation} of a distribution (a sample) with respect to a point in which depth is calculated. In their approach instead of a distribution  \eqn{ {P}^{X} } , a distribution  \eqn{ {{P}_{x}}=1/2{{P}^{X}}+1/2{{P}^{2x-X}} }  is used. For any  \eqn{ \beta \in [0,1] } , let us introduce the smallest depth region bigger or equal to  \eqn{ \beta  } ,  \deqn{ {R}^{\beta }(F)=\bigcap\limits_{\alpha \in A(\beta )}{{{D}_{\alpha }}}(F), }  where  \eqn{ A(\beta )=\left\{ \alpha \ge 0:P\left[ {{D}_{\alpha }}(F) \right]\ge \beta  \right\} } . Then for a locality parameter  \eqn{ \beta  }  we can take a neighbourhood of a point  \eqn{ x }  as  \eqn{ R_{x}^{\beta }(P) } .
#'  
#'  Formally, let  \eqn{ D(\cdot,P) }  be a depth function. Then the {local depth }with the locality parameter  \eqn{ \beta  }  and w.r.t. a point  \eqn{ x }  is defined as  \deqn{ L{{D}^{\beta }}(z,P):z\to D(z,P_{x}^{\beta }), }  where  \eqn{ P_{x}^{\beta }(\cdot )=P\left( \cdot |R_{x}^{\beta }(P) \right) }  is cond. distr. of  \eqn{ P }  conditioned on  \eqn{ R_{x}^{\beta }(P) } .
#'  
#'  @references 
#'  
#'  Paindaveine, D., Van Bever, G. (2013) From depth to local depth : a focus on centrality. Journal of the American Statistical Association 105, 1105-1119 (2013).
#'  
#'  @examples
#' \dontrun{
#' # EXAMPLE 1
#' data = mvrnorm(100, c(0,5), diag(2)*5)
#' #by default depth2 = depth1
#' depthLocal(data, data, depth1 = "LP")
#' depthLocal(data, data, depth1 = "LP", depth2 = "Projection")
#' ## Depthcontour
#' depthContour(data, method = "Local", depth1 = "LP")
#' # EXAMPLE 2
#' data(inf.mort,maesles.imm)
#' data1990=na.omit(cbind(inf.mort[,1],maesles.imm[,1]))
#' depthContour(data1990, method = "Local", depth1 = "LP",beta=0.3)
#' 
#' #EXAMPLE 3
#' Sigma1 = matrix(c(10,3,3,2),2,2)
#' X1 = mvrnorm(n= 8500, mu= c(0,0),Sigma1)
#' Sigma2 = matrix(c(10,0,0,2),2,2)
#' X2 = mvrnorm(n= 1500, mu= c(-10,6),Sigma2)
#' BALLOT=rbind(X1,X2)
#' 
#' train <- sample(1:10000, 100)
#' data<-BALLOT[train,]
#' depthContour(data, method = "Local", depth1 = "Projection",beta=0.3)
#' }
#' 
depthLocal = function(u, X, beta=0.5,
                      depth1="Projection", depth2=depth1, name = "X", ...) 
{
  if(missing(X)) X = u
  
  depths = 1:nrow(u)
  for(i in 1:nrow(u)) depths[i] = .depthLocal(u[i,,drop = FALSE], X, beta, depth1, depth2, ...)
  
  new("DepthLocal", depths, u = u, X = X, method = "Local", name = name,depth1 = depth1, depth2 = depth2)
  
}
 

