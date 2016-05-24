#'@title SIMPLS-CA: SIMPLS Canonical Analysis
#'
#'@description
#'The function \code{simplsca} performs the SIMPLS Canonical Analysis 
#'algorithm as described in Michel Tenenhaus book \emph{La Regression PLS}, chapter 5.
#'
#'@details
#'No missing data are allowed.
#'
#'@param X Numeric matrix or data frame with two or more columns (X-block).
#'@param Y Numeric matrix or data frame with two or more columns (Y-block).
#'@param comps Number of components to be extracted.
#'(\code{TRUE} by default).
#'@return An object of class \code{"simplsca"}, basically a list with the
#'following elements:
#'@return \item{x.scores}{scores of the X-block (also known as T
#'components)}
#'@return \item{x.wgs}{weights of the X-block}
#'@return \item{y.scores}{scores of the Y-block (also known as U
#'components)}
#'@return \item{y.wgs}{weights of the Y-block}
#'@return \item{cor.xt}{correlations between X and T}
#'@return \item{cor.yu}{correlations between Y and U}
#'@return \item{cor.xu}{correlations between X and U}
#'@return \item{cor.yt}{correlations between Y and T}
#'@return \item{cor.tu}{correlations between T and U}
#'@return \item{R2XT}{explained variance of X by T}
#'@return \item{R2YT}{explained variance of Y by T}
#'@return \item{R2YU}{explained variance of Y by U}
#'@return \item{R2XU}{explained variance of X by U}
#'@author Gaston Sanchez
#'@seealso \code{\link{plot.simplsca}}, \code{\link{simpls}}
#'@references Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et
#'Pratique.} Paris: Editions TECHNIP.
#'@export
#'@examples
#'
#'  \dontrun{
#'  # load data linnerud
#'  data(linnerud)
#'
#'  # apply inter-battery method
#'  my_simca = simplsca(linnerud[,1:3], linnerud[,4:6])
#'
#'  # plot variables
#'  plot(my_simca, what="variables")
#'  }
#'
simplsca <- function(X, Y, comps=2)
{
  # =======================================================
  # checking arguments
  # =======================================================
  # matrix X
  if (is.null(dim(X))) stop("\nX is not a matrix")
  if (!is.matrix(X)) X = as.matrix(X)
  p = ncol(X)
  if (p == 1) stop("\nX must have more than one column")
  if (any(!is.finite(X)))
    stop("\ninfinite, NA or NaN values in X")
  # matrix Y  
  if (is.null(dim(Y))) stop("\nY is not a matrix")
  if (!is.matrix(Y)) Y = as.matrix(Y)
  q = ncol(Y)
  if (q == 1) stop("\nY must have more than one column")
  if (any(!is.finite(Y)))
    stop("\ninfinite, NA or NaN values in Y") 
  n = nrow(X)
  if (nrow(Y) != n)
    stop("\nX and Y have different number of rows")
  # names
  if (is.null(colnames(X)))
    colnames(X) = paste(rep("X",p), 1:p, sep="")
  if (is.null(rownames(X)))
    rownames(X) = rep(1:n)
  if (is.null(colnames(Y)))
    colnames(Y) = paste(rep("Y",q), 1:q, sep="")
  if (is.null(rownames(Y)))
    rownames(Y) = rep(1:n)        
  # number of components 
  if (mode(comps)!="numeric" || length(comps)!=1 || 
    comps<=1 || (comps%%1)!=0)
    comps = 2
  svdX = svd(X, nu = 0)
  Xrank = sum(svdX$d > 0.0001^2)
  if (Xrank == 0) stop("\nrank of X = 0: variables are numerically const")
  svdY = svd(Y, nu = 0)
  Yrank = sum(svdY$d > 0.0001^2)
  if (Yrank == 0) stop("\nrank of Y = 0: variables are numerically const")
  nc = min(comps, min(Xrank, Yrank))
  # centering and normalizing
  Xold = scale(X)
  Yold = scale(Y)
  
  # =======================================================
  # SIMPLS-CA algorithm
  # =======================================================
  A = matrix(0, p, nc)
  B = matrix(0, q, nc)
  Xt = matrix(0, n, nc)
  Yu = matrix(0, n, nc)
  for (h in 1:nc)
  {
    # singular value decomposition
    XY_svd = svd(t(Xold) %*% Yold)
    a.new = XY_svd$u[,1]
    b.new = XY_svd$v[,1]
    A[,h] = a.new
    B[,h] = b.new
    # components t's and u's
    t.new = Xold %*% a.new
    u.new = Yold %*% b.new
    Xt[,h] = t.new
    Yu[,h] = u.new
    # Orthogonal projectors
    cx = t(Xold) %*% t.new
    Px = cx %*% solve(t(cx) %*% cx) %*% t(cx)
    cy = t(Yold) %*% u.new
    Py = cy %*% solve(t(cy) %*% cy) %*% t(cy)
    # deflate
    Xold = Xold - (Xold %*% Px)
    Yold = Yold - (Yold %*% Py)
  }  
  dimnames(A) = list(colnames(X), paste(rep("t",nc),1:ncol(A),sep=""))
  dimnames(B) = list(colnames(Y), paste(rep("u",nc),1:ncol(B),sep=""))
  dimnames(Xt) = list(rownames(X), paste(rep("t",nc),1:ncol(A),sep=""))
  dimnames(Yu) = list(rownames(Y), paste(rep("u",nc),1:ncol(B),sep=""))
  
  # =======================================================
  # Complementary results
  # =======================================================
  # correlations between variables and components
  cor.xt = cor(X, Xt)
  cor.xt = cor(X, Xt)
  cor.yu = cor(Y, Yu)
  cor.xu = cor(X, Yu)
  cor.yt = cor(Y, Xt)
  cor.tu = cor(cbind(Xt, Yu))
  # explained variance
  R2XT = t(apply(cor.xt^2, 1, cumsum))
  R2YT = t(apply(cor.yt^2, 1, cumsum))
  R2YU = t(apply(cor.yu^2, 1, cumsum))
  R2XU = t(apply(cor.xu^2, 1, cumsum))
  
  # results 
  res = list(x.scores = Xt,
             x.wgs = A,
             y.scores = Yu,
             y.wgs = B, 
             cor.xt = cor.xt, 
             cor.yu = cor.yu, 
             cor.xu = cor.xu, 
             cor.yt = cor.yt,
             cor.tu = cor.tu,
             R2XT = R2XT,
             R2YT = R2YT,
             R2YU = R2YU,
             R2XU = R2XU)
  class(res) = "simplsca"
  return(res)
}
