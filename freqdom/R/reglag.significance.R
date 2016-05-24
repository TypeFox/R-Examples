#' Check significance of estimated coefficents in \code{\link{speclagreg}} estimator
#'
#' @title Test significance of coefficients in linear model estimator
#' @param X first multivariate time series
#' @param Y second multivariate time series
#' @param A estimated operators
#' @param alpha significance level
#' @param plot plot the results
#' @param ... arguments passed to \code{\link{plot}} function
#' @return list with a quantile and test statistics for each lag
#' @examples
#' n = 200
#' d = 5
#' X = rar(n,d=d,Psi=matrix(0,d,d))  			# independent d-dim variables
#' w = 0.4
#' Y = w*X + (1-w)*rar(n,d=d,Psi=matrix(0,d,d))	# independent d-dim variables
#' A = speclagreg(X, Y, lags=-2:2)
#' W = reglag.significance(X, Y, A, alpha = 0.05)
#' @importFrom graphics abline plot title
#' @importFrom stats pchisq qchisq quantile
#' @export 
reglag.significance = function(X, Y, A, alpha = 0.05, plot = FALSE, ...)
{
  if (!is.matrix(Y))
    stop("Y must be a matrix")
  if (!is.matrix(X))
    stop("X must be a matrix")
  if (!is.timedom(A))
    stop("A must be a time domain operator")
  if (alpha > 1 || alpha < 0)
    stop("alpha must be in [0,1]")
  
  dX = dim(X)[2]
  dY = dim(Y)[2]
  n = dim(X)[1]
  lags = A$lags
  
  # Compute spectral densities
  SX = spectral.density(X)
  residua = Y - A %c% X
  
  SY = spectral.density(Y)
  C = lagged.cov(X,lag=0)
  
  Cinv = pseudoinverse(C,K = dim(C)[1])
  #Cinv = solve(C)
  
  for (i in 1:length(SX$freq)){
    SX$operators[i,,] = Cinv %*% SX$operators[i,,] %*% Cinv
  }

  brillinger = TRUE
  
  if (brillinger){
    # Get the asymptotic distribution of operators under the null
    # as Brillinger suggests in Theorem 8.10.2
    PROD = freqdom.kronecker(SX,SY)
    B = invfourier(PROD)$operators[1,,] / n
    
    invB = pseudoinverse(B,K = dim(B)[1])
    
    # 'normalise' operators with invB
    W = c()
    for (i in 1:length(A$lags))
      W = c(W,c(A$operators[i,,]) %*% invB %*% c(A$operators[i,,]))
    W = Re(W)
  }
  else{ # TODO: NOT READY YET
    # Make it a bit simpler
    n = dim(X)[1]
    PROD = freqdom.ratio(SY,SX,n)
    B = invfourier(PROD)$operators[1,,] / n
    invB = solve(B) # Invert B
    
    # 'normalise' operators with invB
    W = c()
    for (i in 1:length(A$lags)){
      res = t(A$operators[i,,]) %*% B %*% (A$operators[i,,])
      W = c(W,trace(res))
    }
  }
  
  
  # the quantile
  q.joint = qchisq(1-alpha/length(A$lags),(2*length(A$lags) - 1) * dX * dY)
  p.joint = pchisq(sum(W),(2*length(A$lags) - 1) * dX * dY)
  cat("Testing lags jointly:\n")
  print(lags)
  cat("stat:",sum(W),"\n")
  cat("p-value:",1-p.joint,"\n")
  
  q = qchisq(1-alpha,dX * dY)
  p = pchisq(W,dX * dX)
  
  res = list()
  res$lags = lags
  res$stat = abs(W)
  res$q.sep = q
  res$p.sep = 1-p
  
  res$q.joint = q.joint
  res$p.joint = 1-p.joint
  
  # Plot the results
  if (plot){
    ylim = c(min(res$stat,res$q.sep,0),max(res$stat,res$q.sep))
  
    args = list(...)
    if (!("ylim" %in% names(args)))
      args$ylim = ylim

    args = c(list(x=res$lags,y=res$stat,ylab="Test statistic value",xlab="lag"),args)

    args$plot = NULL
    args$alpha = NULL
    
    do.call(graphics::plot,args)
    abline(h=res$q.sep,col=3)
    title("Significance of lags (each one separetly)")
  }
  
  res
}
