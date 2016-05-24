#'@title Tucker's Inter-Battery Method of Factor Analysis
#'
#'@description
#'The function \code{interbat} performs Tucker's Inter-Battery method of factor analysis as described 
#'in Michel Tenenhaus book \emph{La Regression PLS}, chapter 3
#'
#'@details
#'Arguments \code{X} and \code{Y} must contain more than one variable.
#'No missing data are allowed.
#'
#'@param X Numeric matrix or data frame with two or more columns (X-block).
#'@param Y Numeric matrix or data frame with two or more columns (Y-block).
#'@param scaled Logical value indicating whether to scale the data
#'(\code{TRUE} by default).
#'@return An object of class \code{"interbat"}, basically a list with the
#'following elements:
#'@return \item{values}{The extracted eigenvalues}
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
#'@return \item{R2X}{explained variance of X by T}
#'@return \item{R2Y}{explained variance of Y by U}
#'@return \item{com.xu}{communality of X with U}
#'@return \item{com.yt}{communality of Y with T}
#'@return \item{statistic}{Phi statistic values for assessing the number 
#'of relevant components}
#'@author Gaston Sanchez
#'@seealso \code{\link{plot.interbat}}, \code{\link{plsca}}
#'@references Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et
#'Pratique.} Paris: Editions TECHNIP.
#'
#'Tucker, L.R. (1958) An inter-battery method of factor analysis. 
#'\emph{Psychometrika}, 23(2): 111-136.
#'@export
#'@examples
#'
#'  \dontrun{
#'  # load data linnerud
#'  data(linnerud)
#'
#'  # apply inter-battery method
#'  ib = interbat(linnerud[,1:3], linnerud[,4:6])
#'
#'  # plot variables
#'  plot(ib, what="variables")
#'
#'  # plot observations
#'  plot(ib, what="observations", comps=c(1,1), where=c("t","u"))
#'  }
#'
interbat <- function(X, Y, scaled=TRUE)
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
  # scaling data
  if (scaled) {
    X = scale(X)
    Y = scale(Y)
  }
  
  # =======================================================
  # Inter-battery method
  # =======================================================
  # correlation matrix
  Rxy = cor(X, Y)
  Ryx = t(Rxy)
  # get eigenvectors
  Rxy_eig = eigen(Rxy %*% Ryx)
  A = Rxy_eig$vectors
  eigs = Rxy_eig$values
  B = eigen(Ryx %*% Rxy)$vectors
  dimnames(A) = list(colnames(X), paste(rep("t",p),1:ncol(A),sep=""))
  dimnames(B) = list(colnames(Y), paste(rep("u",q),1:ncol(B),sep=""))
  # get components t's and u's
  Xt = X %*% A
  Yu = Y %*% B
  dimnames(Xt) = list(rownames(X), paste(rep("t",p),1:ncol(A),sep=""))
  dimnames(Yu) = list(rownames(Y), paste(rep("u",q),1:ncol(B),sep=""))

  # =======================================================
  # Complementary results
  # =======================================================
  # correlations between variables and scores
  cor.xt = cor(X, Xt)
  cor.yu = cor(Y, Yu)
  cor.xu = cor(X, Yu)
  cor.yt = cor(Y, Xt)
  cor.tu = cor(cbind(Xt, Yu))
  # explained variance
  R2X = t(apply(cor.xt^2, 1, cumsum))
  R2Y = t(apply(cor.yu^2, 1, cumsum))
  # communalities
  Com.XU = t(apply(cor.xu^2, 1, cumsum))
  Com.YT = t(apply(cor.yt^2, 1, cumsum))
  # statistic for assessing the number of relevant components
  m = min(p, q)
  varxt = apply(Xt, 2, var)
  varyu = apply(Yu, 2, var)
  stat = rep(0, m)  # statistic vector
  dfs = rep(0, m)  # degrees of freedom
  stat[1] = n * sum(eigs)
  dfs[1] = p * q
  for (i in 1:(m-1))
  {
    k = n * (p - i) * (q - i)
    num = k * sum(eigs[-(1:i)])
    denum = (p - sum(varxt[1:i])) * (q - sum(varyu[1:i]))
    stat[i+1] = num / denum
    dfs[i+1] = (p - i) * (q - i)
  }
  pvals = 1 - pchisq(stat, dfs)
  stat_df = data.frame(phi=stat, df=dfs, p.value=pvals)
  # results 
  res = list(values = eigs, 
             x.scores = Xt, 
             x.wgs = A,
             y.scores = Yu,
             y.wgs = B,
             cor.xt = cor.xt, 
             cor.yu = cor.yu, 
             cor.xu = cor.xu, 
             cor.yt = cor.yt,
             cor.tu = cor.tu, 
             R2X = R2X, 
             R2Y = R2Y, 
             com.xu = Com.XU, 
             com.yt = Com.YT,
             statistic = stat_df)
  class(res) = "interbat"
  return(res)
}
